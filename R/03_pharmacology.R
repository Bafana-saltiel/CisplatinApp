# ============================================================
# 4. CORE PHARMACOLOGICAL FUNCTIONS
# ============================================================

make_seed <- function(...) {
  seed_string <- paste(..., sep = "_")
  ints <- utf8ToInt(seed_string)
  seed_value <- sum(ints * seq_along(ints))
  abs(seed_value %% .Machine$integer.max)
}

predict_viability <- function(dose, ec50, lower, upper, slope) {
  lower + (upper - lower) / (1 + (dose / ec50)^slope)
}

predict_inhibition <- function(dose, ec50, lower, upper, slope) {
  viability <- predict_viability(dose, ec50, lower, upper, slope)
  inhibition <- (1 - viability) * 100
  pmax(pmin(inhibition, 100), 0)
}

bliss_multi <- function(E) {
  E <- pmax(pmin(as.numeric(E), 100), 0) / 100
  100 * (1 - prod(1 - E))
}

hsa_multi <- function(E) {
  E <- pmax(pmin(as.numeric(E), 100), 0)
  max(E, na.rm = TRUE)
}

# Invert one fitted monotherapy curve to estimate the dose required to
# produce a requested inhibition level. This supports a generalized
# Loewe concentration-additivity calculation.
inverse_inhibition_dose <- function(
    target_inhibition,
    ec50,
    lower,
    upper,
    slope
) {
  pars <- c(target_inhibition, ec50, lower, upper, slope)

  if (any(!is.finite(pars)) || ec50 <= 0 || slope == 0) {
    return(NA_real_)
  }

  target_inhibition <- pmax(pmin(target_inhibition, 99.999), 0.001)
  target_viability <- 1 - target_inhibition / 100

  denominator <- target_viability - lower
  numerator <- upper - target_viability

  if (
    denominator <= 0 ||
    numerator <= 0 ||
    !is.finite(numerator / denominator)
  ) {
    return(NA_real_)
  }

  dose <- ec50 * (numerator / denominator)^(1 / slope)

  if (!is.finite(dose) || dose <= 0) {
    NA_real_
  } else {
    dose
  }
}

loewe_expected_multi <- function(doses, drug_parameters) {
  doses <- as.numeric(doses)

  if (
    length(doses) == 0 ||
    nrow(drug_parameters) != length(doses) ||
    any(!is.finite(doses))
  ) {
    return(NA_real_)
  }

  loewe_equation <- function(effect) {
    equivalent_doses <- vapply(
      seq_len(nrow(drug_parameters)),
      function(i) {
        inverse_inhibition_dose(
          target_inhibition = effect,
          ec50 = drug_parameters$EC50[[i]],
          lower = drug_parameters$LowerAsymptote[[i]],
          upper = drug_parameters$UpperAsymptote[[i]],
          slope = drug_parameters$Slope[[i]]
        )
      },
      numeric(1)
    )

    valid <- is.finite(equivalent_doses) & equivalent_doses > 0

    if (!any(valid)) {
      return(NA_real_)
    }

    sum(doses[valid] / equivalent_doses[valid]) - 1
  }

  lower_effect <- 0.001
  upper_effect <- 99.999

  lower_value <- suppressWarnings(loewe_equation(lower_effect))
  upper_value <- suppressWarnings(loewe_equation(upper_effect))

  if (
    !is.finite(lower_value) ||
    !is.finite(upper_value) ||
    lower_value * upper_value > 0
  ) {
    individual_effects <- vapply(
      seq_len(nrow(drug_parameters)),
      function(i) {
        predict_inhibition(
          dose = doses[[i]],
          ec50 = drug_parameters$EC50[[i]],
          lower = drug_parameters$LowerAsymptote[[i]],
          upper = drug_parameters$UpperAsymptote[[i]],
          slope = drug_parameters$Slope[[i]]
        )
      },
      numeric(1)
    )

    return(hsa_multi(individual_effects))
  }

  root <- tryCatch(
    uniroot(
      loewe_equation,
      interval = c(lower_effect, upper_effect),
      tol = 1e-6
    )$root,
    error = function(e) NA_real_
  )

  pmax(pmin(root, 100), 0)
}

# ZIP normally estimates zero-interaction potency from a complete response
# surface. The app reconstructs responses from single-agent fitted curves,
# so this implementation uses a transparent curve-based ZIP approximation:
# the midpoint of Bliss independence and Loewe concentration additivity.
zip_expected_multi <- function(E, doses, drug_parameters) {
  bliss_value <- bliss_multi(E)
  loewe_value <- loewe_expected_multi(doses, drug_parameters)

  if (!is.finite(loewe_value)) {
    return(bliss_value)
  }

  mean(c(bliss_value, loewe_value))
}

calculate_expected_response <- function(
    inhibition,
    doses,
    drug_parameters,
    method
) {
  method <- match.arg(
    method,
    c("bliss", "hsa", "loewe", "zip")
  )

  switch(
    method,
    bliss = bliss_multi(inhibition),
    hsa = hsa_multi(inhibition),
    loewe = loewe_expected_multi(doses, drug_parameters),
    zip = zip_expected_multi(
      E = inhibition,
      doses = doses,
      drug_parameters = drug_parameters
    )
  )
}

synergy_method_label <- function(method) {
  switch(
    method,
    bliss = "Bliss independence",
    hsa = "Highest Single Agent (HSA)",
    loewe = "Loewe additivity",
    zip = "ZIP approximation",
    method
  )
}

safe_mean <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_min <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)

trapz_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

find_drug_rows <- function(drug_name, df = Cervix_GDSC_All) {
  drug_name_upper <- stringr::str_to_upper(as.character(drug_name)[1])

  compound_id <- switch(
    drug_name_upper,
    "CISPLATIN" = "DPC-001793",
    "C225" = "DPC-001418",
    "CETUXIMAB" = "DPC-001418",
    "PACLITAXEL" = "DPC-004880",
    "TOPOTECAN" = "DPC-006564",
    "FLUOROURACIL" = "DPC-002828",
    "METHOTREXATE" = "DPC-000234",
    "BLEOMYCIN" = "DPC-001132",
    "VINORELBINE" = "DPC-002684",
    "VINBLASTINE" = "DPC-006819",
    "VELIPARIB" = "DPC-000236",
    NA_character_
  )

  if (is.na(compound_id)) return(tibble::tibble())

  df %>%
    left_join(
      PortalCompounds %>%
        select(CompoundID, Synonyms, ChEMBLID, PubChemCID, TargetOrMechanism),
      by = "CompoundID"
    ) %>%
    filter(CompoundID == compound_id)
}

simulate_reference_combo <- function(model_id, combo_name, therapeutic_class, drug_names, doses, noise_sd, synergy_method) {
  df_model <- Cervix_GDSC_All %>% filter(ModelID == model_id)

  matched_drugs <- map_dfr(drug_names, function(drug) {
    rows <- find_drug_rows(drug, df_model)
    if (nrow(rows) == 0) return(tibble())
    rows %>% slice(1) %>% mutate(QueryDrug = drug)
  })

  if (nrow(matched_drugs) != length(drug_names)) return(tibble())

  dose_grid <- expand.grid(replicate(length(drug_names), doses, simplify = FALSE))
  colnames(dose_grid) <- paste0(matched_drugs$QueryDrug, "_Dose")

  inhibition_matrix <- map_dfc(seq_len(nrow(matched_drugs)), function(i) {
    inhib <- predict_inhibition(
      dose = dose_grid[[i]],
      ec50 = matched_drugs$EC50[[i]],
      lower = matched_drugs$LowerAsymptote[[i]],
      upper = matched_drugs$UpperAsymptote[[i]],
      slope = matched_drugs$Slope[[i]]
    )
    tibble(!!paste0(matched_drugs$QueryDrug[[i]], "_Response") := inhib)
  })

  expected_response <- vapply(
    seq_len(nrow(dose_grid)),
    function(i) {
      calculate_expected_response(
        inhibition = as.numeric(inhibition_matrix[i, , drop = TRUE]),
        doses = as.numeric(dose_grid[i, , drop = TRUE]),
        drug_parameters = matched_drugs,
        method = synergy_method
      )
    },
    numeric(1)
  )

  # A common observed-response surface is generated independently of the
  # selected scoring model so that users can compare models consistently.
  bliss_baseline <- apply(inhibition_matrix, 1, bliss_multi)

  set.seed(make_seed(model_id, combo_name, "REFERENCE"))
  E_observed <- pmax(
    pmin(
      bliss_baseline + rnorm(length(bliss_baseline), mean = 0, sd = noise_sd),
      100
    ),
    0
  )

  bind_cols(
    tibble(
      ModelID = model_id,
      CellLineName = matched_drugs$CellLineName[[1]],
      DepmapModelType = matched_drugs$DepmapModelType[[1]],
      ClinicalCombination = combo_name,
      TherapeuticClass = therapeutic_class,
      RequiredDrugs = paste(drug_names, collapse = " + "),
      MatchedCompoundNames = paste(matched_drugs$CompoundName, collapse = " + "),
      MatchedCompoundIDs = paste(matched_drugs$CompoundID, collapse = " + ")
    ),
    dose_grid,
    inhibition_matrix,
    tibble(
      SynergyMethod = synergy_method_label(synergy_method),
      ExpectedResponse = expected_response,
      Observed = E_observed,
      SynergyScore = E_observed - expected_response
    )
  )
}

evaluate_monotherapy <- function(model_id, partner, doses) {
  if (is.null(partner) || nrow(partner) == 0) return(tibble())

  row <- partner %>% slice(1)
  required_values <- c(row$EC50[[1]], row$LowerAsymptote[[1]], row$UpperAsymptote[[1]], row$Slope[[1]])
  if (any(is.na(required_values))) return(tibble())

  dose_values <- sort(unique(as.numeric(doses)))
  dose_values <- dose_values[is.finite(dose_values)]
  if (length(dose_values) == 0) return(tibble())

  inhibition <- predict_inhibition(
    dose = dose_values,
    ec50 = row$EC50[[1]],
    lower = row$LowerAsymptote[[1]],
    upper = row$UpperAsymptote[[1]],
    slope = row$Slope[[1]]
  )

  log_dose <- log10(dose_values)
  auc_raw <- trapz_auc(log_dose, inhibition)
  dose_range <- diff(range(log_dose))
  auc_normalised <- if (is.finite(auc_raw) && is.finite(dose_range) && dose_range > 0) {
    auc_raw / dose_range
  } else {
    safe_mean(inhibition)
  }

  mean_response <- safe_mean(inhibition)
  monotherapy_class <- case_when(
    is.na(mean_response) ~ "Not available",
    mean_response >= 75 ~ "High intrinsic response",
    mean_response >= 40 ~ "Moderate intrinsic response",
    mean_response > 0 ~ "Low intrinsic response",
    TRUE ~ "No predicted response"
  )

  tibble(
    ModelID = model_id,
    CellLineName = row$CellLineName[[1]],
    DepmapModelType = row$DepmapModelType[[1]],
    ThirdDrug = row$CompoundName[[1]],
    ThirdDrugCompoundID = row$CompoundID[[1]],
    Dose = dose_values,
    PredictedMonotherapyInhibition = inhibition,
    EC50 = row$EC50[[1]],
    LowerAsymptote = row$LowerAsymptote[[1]],
    UpperAsymptote = row$UpperAsymptote[[1]],
    Slope = row$Slope[[1]],
    MinimumDose = row$MinimumDose[[1]],
    MaximumDose = row$MaximumDose[[1]],
    DoseUnit = row$DoseUnit[[1]],
    MeanMonotherapyResponse = mean_response,
    MaxMonotherapyResponse = safe_max(inhibition),
    MinMonotherapyResponse = safe_min(inhibition),
    SDMonotherapyResponse = safe_sd(inhibition),
    MonotherapyAUC = auc_normalised,
    MonotherapyClass = monotherapy_class
  )
}

simulate_reference_plus_partner <- function(model_id, combo_name, therapeutic_class, drug_names, partner, doses, noise_sd, synergy_method) {
  df_model <- Cervix_GDSC_All %>% filter(ModelID == model_id)

  matched_reference <- map_dfr(drug_names, function(drug) {
    rows <- find_drug_rows(drug, df_model)
    if (nrow(rows) == 0) return(tibble())
    rows %>% slice(1) %>% mutate(QueryDrug = drug)
  })

  if (nrow(matched_reference) != length(drug_names)) return(tibble())

  matched_drugs <- bind_rows(
    matched_reference,
    partner %>% mutate(QueryDrug = paste0("ThirdDrug_", CompoundName))
  )

  dose_grid <- expand.grid(replicate(nrow(matched_drugs), doses, simplify = FALSE))
  colnames(dose_grid) <- paste0(matched_drugs$QueryDrug, "_Dose")

  inhibition_matrix <- map_dfc(seq_len(nrow(matched_drugs)), function(i) {
    inhib <- predict_inhibition(
      dose = dose_grid[[i]],
      ec50 = matched_drugs$EC50[[i]],
      lower = matched_drugs$LowerAsymptote[[i]],
      upper = matched_drugs$UpperAsymptote[[i]],
      slope = matched_drugs$Slope[[i]]
    )
    tibble(!!paste0(matched_drugs$QueryDrug[[i]], "_Response") := inhib)
  })

  expected_response <- vapply(
    seq_len(nrow(dose_grid)),
    function(i) {
      calculate_expected_response(
        inhibition = as.numeric(inhibition_matrix[i, , drop = TRUE]),
        doses = as.numeric(dose_grid[i, , drop = TRUE]),
        drug_parameters = matched_drugs,
        method = synergy_method
      )
    },
    numeric(1)
  )

  bliss_baseline <- apply(inhibition_matrix, 1, bliss_multi)

  set.seed(make_seed(model_id, combo_name, partner$CompoundID[[1]]))
  E_observed <- pmax(
    pmin(
      bliss_baseline + rnorm(length(bliss_baseline), mean = 0, sd = noise_sd),
      100
    ),
    0
  )

  bind_cols(
    tibble(
      ModelID = model_id,
      CellLineName = matched_reference$CellLineName[[1]],
      DepmapModelType = matched_reference$DepmapModelType[[1]],
      ReferenceCombination = combo_name,
      TherapeuticClass = therapeutic_class,
      ThirdDrug = partner$CompoundName[[1]],
      ThirdDrugCompoundID = partner$CompoundID[[1]],
      ReferenceCompoundNames = paste(matched_reference$CompoundName, collapse = " + "),
      ReferenceCompoundIDs = paste(matched_reference$CompoundID, collapse = " + ")
    ),
    dose_grid,
    inhibition_matrix,
    tibble(
      SynergyMethod = synergy_method_label(synergy_method),
      ExpectedResponse = expected_response,
      Observed = E_observed,
      SynergyScore = E_observed - expected_response
    )
  )
}

run_reference_pipeline <- function(combo_name, drug_names, compound_ids, therapeutic_class, doses, noise_sd, synergy_method, selected_partner_ids = NULL) {
  eligible_models <- Cervix_GDSC_All %>%
    filter(CompoundID %in% compound_ids) %>%
    distinct(ModelID, CellLineName, DepmapModelType, CompoundID) %>%
    count(ModelID, CellLineName, DepmapModelType, name = "n_reference_drugs") %>%
    filter(n_reference_drugs == length(unique(compound_ids))) %>%
    select(ModelID, CellLineName, DepmapModelType)

  empty_result <- function() {
    list(
      eligible_models = eligible_models,
      reference_synergy = tibble(),
      reference_summary = tibble(),
      monotherapy_dose_response = tibble(),
      monotherapy_summary = tibble(),
      third_synergy = tibble(),
      third_ranking = tibble(),
      top10 = tibble()
    )
  }

  if (nrow(eligible_models) == 0) return(empty_result())

  reference_synergy <- map_dfr(
    eligible_models$ModelID,
    ~ simulate_reference_combo(
      model_id = .x,
      combo_name = combo_name,
      therapeutic_class = therapeutic_class,
      drug_names = drug_names,
      doses = doses,
      noise_sd = noise_sd,
      synergy_method = synergy_method
    )
  )

  if (nrow(reference_synergy) == 0) return(empty_result())

  reference_summary <- reference_synergy %>%
    group_by(
      ModelID, CellLineName, DepmapModelType,
      ClinicalCombination, TherapeuticClass, SynergyMethod,
      RequiredDrugs, MatchedCompoundNames, MatchedCompoundIDs
    ) %>%
    summarise(
      ReferenceMeanSynergy = safe_mean(SynergyScore),
      ReferenceMaxSynergy = safe_max(SynergyScore),
      ReferenceMinSynergy = safe_min(SynergyScore),
      ReferenceSDSynergy = safe_sd(SynergyScore),
      ReferenceMeanExpected = safe_mean(ExpectedResponse),
      ReferenceMeanObserved = safe_mean(Observed),
      .groups = "drop"
    )

  run_model <- function(model_id) {
    df_model <- Cervix_GDSC_All %>% filter(ModelID == model_id)

    partners <- df_model %>%
      filter(!CompoundID %in% compound_ids) %>%
      distinct(
        CompoundID, CompoundName, ModelID, CellLineName, DepmapModelType,
        EC50, LowerAsymptote, UpperAsymptote, Slope,
        MinimumDose, MaximumDose, DoseUnit,
        .keep_all = TRUE
      ) %>%
      filter(!is.na(EC50), !is.na(LowerAsymptote), !is.na(UpperAsymptote), !is.na(Slope))

    if (!is.null(selected_partner_ids)) {
      partners <- partners %>% filter(CompoundID %in% selected_partner_ids)
    }

    if (nrow(partners) == 0) return(list(monotherapy = tibble(), third_synergy = tibble()))

    monotherapy <- map_dfr(
      seq_len(nrow(partners)),
      ~ evaluate_monotherapy(model_id, partners[.x, , drop = FALSE], doses)
    )

    third_synergy <- map_dfr(
      seq_len(nrow(partners)),
      ~ simulate_reference_plus_partner(
        model_id = model_id,
        combo_name = combo_name,
        therapeutic_class = therapeutic_class,
        drug_names = drug_names,
        partner = partners[.x, , drop = FALSE],
        doses = doses,
        noise_sd = noise_sd,
        synergy_method = synergy_method
      )
    )

    list(monotherapy = monotherapy, third_synergy = third_synergy)
  }

  model_results <- map(eligible_models$ModelID, run_model)
  monotherapy_dose_response <- map_dfr(model_results, "monotherapy")
  third_synergy <- map_dfr(model_results, "third_synergy")

  monotherapy_summary <- monotherapy_dose_response %>%
    distinct(
      ModelID, CellLineName, DepmapModelType,
      ThirdDrug, ThirdDrugCompoundID,
      EC50, LowerAsymptote, UpperAsymptote, Slope,
      MinimumDose, MaximumDose, DoseUnit,
      MeanMonotherapyResponse, MaxMonotherapyResponse,
      MinMonotherapyResponse, SDMonotherapyResponse,
      MonotherapyAUC, MonotherapyClass
    ) %>%
    arrange(CellLineName, desc(MeanMonotherapyResponse))

  if (nrow(third_synergy) == 0) {
    return(list(
      eligible_models = eligible_models,
      reference_synergy = reference_synergy,
      reference_summary = reference_summary,
      monotherapy_dose_response = monotherapy_dose_response,
      monotherapy_summary = monotherapy_summary,
      third_synergy = tibble(),
      third_ranking = tibble(),
      top10 = tibble()
    ))
  }

  third_ranking <- third_synergy %>%
    group_by(
      ModelID, CellLineName, DepmapModelType,
      ReferenceCombination, TherapeuticClass, SynergyMethod,
      ThirdDrug, ThirdDrugCompoundID,
      ReferenceCompoundNames, ReferenceCompoundIDs
    ) %>%
    summarise(
      n_dose_combinations = n(),
      MeanSynergy = safe_mean(SynergyScore),
      MaxSynergy = safe_max(SynergyScore),
      MinSynergy = safe_min(SynergyScore),
      SDSynergy = safe_sd(SynergyScore),
      MeanExpected = safe_mean(ExpectedResponse),
      MeanObserved = safe_mean(Observed),
      .groups = "drop"
    ) %>%
    left_join(reference_summary %>% select(ModelID, ReferenceMeanSynergy), by = "ModelID") %>%
    left_join(
      monotherapy_summary %>%
        select(
          ModelID, ThirdDrugCompoundID, ThirdDrug, EC50,
          MeanMonotherapyResponse, MaxMonotherapyResponse,
          MinMonotherapyResponse, SDMonotherapyResponse,
          MonotherapyAUC, MonotherapyClass
        ),
      by = c("ModelID", "ThirdDrugCompoundID", "ThirdDrug")
    ) %>%
    mutate(
      Delta_vs_Reference = MeanSynergy - ReferenceMeanSynergy,
      CombinationObservedGain_vs_Monotherapy = MeanObserved - MeanMonotherapyResponse,
      CombinationExpectedGain_vs_Monotherapy = MeanExpected - MeanMonotherapyResponse,
      PharmacologicalInterpretation = case_when(
        MeanSynergy > ReferenceMeanSynergy & MeanMonotherapyResponse < 40 ~ "Combination-specific candidate",
        MeanSynergy > ReferenceMeanSynergy & MeanMonotherapyResponse >= 75 ~ "Strong monotherapy with additional combination benefit",
        MeanSynergy > ReferenceMeanSynergy & MeanMonotherapyResponse >= 40 ~ "Moderate monotherapy with combination benefit",
        MeanSynergy <= ReferenceMeanSynergy & MeanMonotherapyResponse >= 75 ~ "Strong monotherapy but limited combination advantage",
        MeanSynergy < 0 & MeanMonotherapyResponse < 40 ~ "Low monotherapy response and antagonistic combination",
        TRUE ~ "Context-dependent or limited combination benefit"
      ),
      Interaction = case_when(
        MeanSynergy > ReferenceMeanSynergy ~ "Stronger_than_reference",
        MeanSynergy < 0 ~ "Antagonistic",
        TRUE ~ "Weaker_or_Additive"
      )
    ) %>%
    arrange(CellLineName, desc(MeanSynergy))

  top10 <- third_ranking %>%
    group_by(CellLineName) %>%
    slice_max(MeanSynergy, n = 10, with_ties = FALSE) %>%
    ungroup()

  list(
    eligible_models = eligible_models,
    reference_synergy = reference_synergy,
    reference_summary = reference_summary,
    monotherapy_dose_response = monotherapy_dose_response,
    monotherapy_summary = monotherapy_summary,
    third_synergy = third_synergy,
    third_ranking = third_ranking,
    top10 = top10
  )
}

