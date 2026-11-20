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

run_measured_screen <- function(reference_ids, candidate_ids, doses,
                                measurements = Measured_Monotherapy) {
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

  list(reference = select(reference, -ReferenceBlissViability),
       ranking = ranking, coverage = select(by_drug, -DoseValues),
       reference_doses = reference_doses)
}
