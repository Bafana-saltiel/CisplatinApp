# Measurements are stored as unclipped vehicle-relative viability fractions.
# Bound them only when used as probabilities in the Bliss reference.
read_measured_viability <- function(path) {
  measurements <- readr::read_csv(
    path, show_col_types = FALSE,
    col_types = readr::cols(ModelID = "c", CompoundID = "c", SampleID = "c",
                           DoseUM = "d", ViabilityFraction = "d",
                           NPlates = "i", NWells = "i")
  )
  columns <- c("ModelID", "CompoundID", "SampleID", "DoseUM", "ViabilityFraction",
               "NPlates", "NWells")
  if (!all(columns %in% names(measurements)) || nrow(readr::problems(measurements))) {
    stop("Invalid measured monotherapy file: ", path)
  }
  if (anyNA(measurements[columns]) ||
      any(!is.finite(measurements$DoseUM) | measurements$DoseUM <= 0) ||
      any(!is.finite(measurements$ViabilityFraction)) ||
      any(measurements$NPlates <= 0 | measurements$NWells < measurements$NPlates) ||
      anyDuplicated(measurements[c("ModelID", "CompoundID", "DoseUM")])) {
    stop("Invalid or duplicate measured monotherapy records: ", path)
  }
  measurements
}

Measured_Monotherapy <- read_measured_viability(
  file.path(DATA_DIR, "GDSC2MeasuredViability_Cervix.csv")
)
Ovarian_Monotherapy <- read_measured_viability(
  file.path(DATA_DIR, "GDSC2MeasuredViability_Ovary.csv")
)
Monotherapy_Fit_Audit <- readr::read_csv(
  file.path(DATA_DIR, "GDSC2MonotherapyFitAudit_Cervix.csv"), show_col_types = FALSE
)

estimate_at_doses <- function(measurements, doses) {
  x <- log10(measurements$DoseUM)
  y <- measurements$ViabilityFraction
  out <- rep(NA_real_, length(doses))
  exact <- match(doses, measurements$DoseUM)
  out[!is.na(exact)] <- y[exact[!is.na(exact)]]
  if (length(x) > 1L) {
    missing <- is.na(exact)
    out[missing] <- stats::approx(x, y, xout = log10(doses[missing]), rule = 1)$y
  }
  out
}

# Fit increasing four-parameter log-logistic inhibition to measured doses.
fit_monotherapy_null <- function(rows) {
  x <- log(rows$DoseUM)
  y <- pmin(1, pmax(0, 1 - rows$ViabilityFraction))
  if (length(x) < 4L || diff(range(y)) < 1e-6) return(NULL)
  predict_curve <- function(par, log_dose) {
    par[1] + (1 - par[1]) * par[2] * plogis(par[4] * (log_dose - par[3]))
  }
  fit <- optim(c(min(y), .8, median(x), 1),
               function(par) sum((y - predict_curve(par, x))^2),
               method = "L-BFGS-B",
               lower = c(0, 0, min(x) - 5, .01),
               upper = c(1, 1, max(x) + 5, 10))
  par <- fit$par
  if (fit$convergence != 0 || (1 - par[1]) * par[2] < 1e-6) return(NULL)
  list(par = par, doses = rows$DoseUM,
       predict = function(dose) predict_curve(par, log(dose)))
}

loewe_on_grid <- function(curves, grid) {
  out <- rep(NA_real_, nrow(grid))
  if (any(vapply(curves, is.null, logical(1)))) return(out)
  bounds <- vapply(curves, function(curve) range(curve$predict(curve$doses)), numeric(2))
  lo <- max(bounds[1, ])
  hi <- min(bounds[2, ])
  if (lo >= hi) return(out)
  dose_sum <- function(effect, doses) {
    total <- rep(0, nrow(doses))
    for (i in seq_along(curves)) {
      p <- curves[[i]]$par
      fraction <- (effect - p[1]) / ((1 - p[1]) * p[2])
      equivalent_dose <- exp(p[3] + qlogis(fraction) / p[4])
      total <- total + doses[[i]] / equivalent_dose
    }
    total
  }
  covered <- which(dose_sum(lo, grid) >= 1 & dose_sum(hi, grid) <= 1)
  grid <- grid[covered, , drop = FALSE]
  lower <- rep(lo, length(covered))
  upper <- rep(hi, length(covered))
  for (step in seq_len(32)) {
    middle <- (lower + upper) / 2
    above <- dose_sum(middle, grid) > 1
    lower <- ifelse(above, middle, lower)
    upper <- ifelse(above, upper, middle)
  }
  out[covered] <- 100 * (lower + upper) / 2
  out
}

mean_supported <- function(values) {
  if (any(is.finite(values))) mean(values, na.rm = TRUE) else NA_real_
}

run_measured_screen <- function(reference_ids, candidate_ids, doses,
                                measurements = Measured_Monotherapy, include_hsa = FALSE,
                                include_loewe = FALSE, include_zip = FALSE) {
  doses <- sort(unique(as.numeric(doses)))
  stopifnot(length(doses) > 0L, all(is.finite(doses) & doses > 0),
            length(reference_ids) >= 2L, !anyDuplicated(reference_ids))

  candidate_ids <- setdiff(candidate_ids, reference_ids)
  by_drug <- measurements %>%
    filter(CompoundID %in% candidate_ids) %>%
    group_by(ModelID, CompoundID) %>%
    summarise(
      MeasuredDoseCount = n(),
      MinimumMeasuredDoseUM = min(DoseUM),
      MaximumMeasuredDoseUM = max(DoseUM),
      DoseValues = list(estimate_at_doses(pick(DoseUM, ViabilityFraction), doses)),
      .groups = "drop"
    ) %>%
    mutate(GridComplete = map_lgl(DoseValues, ~ all(is.finite(.x))))

  complete <- by_drug %>%
    filter(GridComplete) %>%
    mutate(
      MeanSingleAgentInhibitionPct = map_dbl(DoseValues, ~ 100 * mean(1 - .x)),
      MeanBoundedViability = map_dbl(DoseValues, ~ mean(pmax(0, pmin(1, .x))))
    )

  # The reference uses the original tested doses for each constituent.
  # Mean Bliss viability over independent dose combinations is the product of means.
  reference_doses <- measurements %>% filter(CompoundID %in% reference_ids)
  reference <- reference_doses %>%
    group_by(ModelID, CompoundID) %>%
    summarise(
      MeanBoundedViability = mean(pmax(0, pmin(1, ViabilityFraction))),
      .groups = "drop"
    ) %>%
    group_by(ModelID) %>%
    summarise(
      ReferenceDrugsCovered = n(),
      ReferenceBlissViability = prod(MeanBoundedViability),
      .groups = "drop"
    ) %>%
    filter(ReferenceDrugsCovered == length(reference_ids)) %>%
    mutate(ReferenceBlissExpectedInhibitionPct = 100 * (1 - ReferenceBlissViability))

  ranking <- complete %>%
    filter(CompoundID %in% candidate_ids) %>%
    inner_join(reference, by = "ModelID") %>%
    mutate(
      ExpandedBlissExpectedInhibitionPct = 100 * (1 - ReferenceBlissViability * MeanBoundedViability),
      BlissExpectedIncrementPP = 100 * ReferenceBlissViability * (1 - MeanBoundedViability)
    ) %>%
    select(ModelID, CompoundID, MeanSingleAgentInhibitionPct,
           ReferenceBlissExpectedInhibitionPct, ExpandedBlissExpectedInhibitionPct,
           BlissExpectedIncrementPP) %>%
    arrange(ModelID, desc(BlissExpectedIncrementPP), CompoundID)

  if (include_hsa) {
    # HSA takes the best single-agent inhibition at each dose tuple.
    # In viability units this is the minimum across component drugs.
    hsa_reference <- reference_doses %>%
      semi_join(reference, by = "ModelID") %>%
      group_by(ModelID) %>%
      summarise(
        BestViability = list(Reduce(
          function(left, right) as.vector(outer(left, right, pmin)),
          split(pmax(0, pmin(1, ViabilityFraction)), CompoundID)
        )),
        .groups = "drop"
      ) %>%
      mutate(ReferenceHSAExpectedInhibitionPct = map_dbl(BestViability, ~ 100 * (1 - mean(.x))))

    hsa_candidates <- complete %>%
      select(ModelID, CompoundID, DoseValues) %>%
      inner_join(hsa_reference, by = "ModelID") %>%
      mutate(
        ExpandedHSAExpectedInhibitionPct = map2_dbl(
          BestViability, DoseValues,
          ~ 100 * (1 - mean(outer(.x, pmax(0, pmin(1, .y)), pmin)))
        ),
        HSAExpectedIncrementPP = ExpandedHSAExpectedInhibitionPct - ReferenceHSAExpectedInhibitionPct
      ) %>%
      select(-BestViability, -DoseValues)

    reference <- left_join(reference, select(hsa_reference, -BestViability), by = "ModelID")
    ranking <- left_join(ranking, hsa_candidates, by = c("ModelID", "CompoundID"))
  }

  if (include_zip) {
    reference$ReferenceZIPExpectedInhibitionPct <- rep(NA_real_, nrow(reference))
    ranking$ExpandedZIPExpectedInhibitionPct <- rep(NA_real_, nrow(ranking))
  }
  if (include_loewe) {
    reference$ReferenceLoeweExpectedInhibitionPct <- rep(NA_real_, nrow(reference))
    reference$ReferenceLoeweCoverage <- rep(NA_real_, nrow(reference))
    ranking$ReferenceLoeweExpectedInhibitionPct <- rep(NA_real_, nrow(ranking))
    ranking$ExpandedLoeweExpectedInhibitionPct <- rep(NA_real_, nrow(ranking))
    ranking$LoeweExpectedIncrementPP <- rep(NA_real_, nrow(ranking))
    ranking$LoeweMatchedCoverage <- rep(NA_real_, nrow(ranking))
  }
  if (include_loewe || include_zip) {
    for (model_id in reference$ModelID) {
      model_rows <- measurements %>% filter(ModelID == model_id)
      model_rank <- which(ranking$ModelID == model_id)
      ids <- c(reference_ids, ranking$CompoundID[model_rank])
      fits <- lapply(split(model_rows, model_rows$CompoundID)[ids], fit_monotherapy_null)
      reference_fits <- fits[seq_along(reference_ids)]
      reference_grid <- lapply(reference_fits, function(curve) curve$doses)
      reference_index <- which(reference$ModelID == model_id)
      if (include_zip && !any(vapply(reference_fits, is.null, logical(1)))) {
        viability <- prod(vapply(reference_fits, function(curve) {
          mean(1 - curve$predict(curve$doses))
        }, numeric(1)))
        reference$ReferenceZIPExpectedInhibitionPct[reference_index] <- 100 * (1 - viability)
        for (j in model_rank) {
          curve <- fits[[ranking$CompoundID[j]]]
          if (!is.null(curve)) {
            expanded <- 100 * (1 - viability * mean(1 - curve$predict(doses)))
            ranking$ExpandedZIPExpectedInhibitionPct[j] <- expanded
          }
        }
      }
      if (include_loewe && !any(vapply(reference_fits, is.null, logical(1)))) {
        grid <- do.call(expand.grid, reference_grid)
        reference_values <- loewe_on_grid(reference_fits, grid)
        reference$ReferenceLoeweExpectedInhibitionPct[reference_index] <- mean_supported(reference_values)
        reference$ReferenceLoeweCoverage[reference_index] <- mean(is.finite(reference_values))
        expanded_grid <- do.call(expand.grid, c(reference_grid, list(doses)))
        matched_reference <- rep(reference_values, length(doses))
        for (j in model_rank) {
          expanded <- loewe_on_grid(c(reference_fits, fits[ranking$CompoundID[j]]), expanded_grid)
          matched <- is.finite(expanded) & is.finite(matched_reference)
          ranking$ReferenceLoeweExpectedInhibitionPct[j] <- mean_supported(matched_reference[matched])
          ranking$ExpandedLoeweExpectedInhibitionPct[j] <- mean_supported(expanded[matched])
          ranking$LoeweExpectedIncrementPP[j] <- mean_supported((expanded - matched_reference)[matched])
          ranking$LoeweMatchedCoverage[j] <- mean(matched)
        }
      }
    }
    if (include_zip) {
      ranking <- ranking %>%
        left_join(select(reference, ModelID, ReferenceZIPExpectedInhibitionPct), by = "ModelID") %>%
        mutate(ZIPExpectedIncrementPP = ExpandedZIPExpectedInhibitionPct - ReferenceZIPExpectedInhibitionPct)
    }
  }

  list(reference = select(reference, -ReferenceBlissViability),
       ranking = ranking, coverage = select(by_drug, -DoseValues),
       reference_doses = reference_doses)
}
