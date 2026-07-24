# ============================================================
# Shiny App: Cisplatin Clinical Reference Combination Explorer
# Updated with monotherapy response, selectable synergy models, anti-angiogenic screening, and validation layers
# ============================================================

required_packages <- c(
  "shiny", "shinythemes", "DT", "dplyr", "stringr",
  "purrr", "tidyr", "ggplot2", "readr", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(DT)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================
# 1. LOAD DATA
# ============================================================

DATA_DIR <- "data"

required_files <- list(
  PortalCompounds = file.path(DATA_DIR, "PortalCompounds.csv"),
  Model = file.path(DATA_DIR, "Model.csv"),
  GDSC2ResponseCurves = file.path(DATA_DIR, "GDSC2ResponseCurves.csv")
)

missing_required <- names(required_files)[!file.exists(unlist(required_files))]

if (length(missing_required) > 0) {
  stop(
    "Missing required data files in data/: ",
    paste(missing_required, collapse = ", ")
  )
}

PortalCompounds <- readr::read_csv(required_files[["PortalCompounds"]], show_col_types = FALSE)
Model <- readr::read_csv(required_files[["Model"]], show_col_types = FALSE)
GDSC2ResponseCurves <- readr::read_csv(required_files[["GDSC2ResponseCurves"]], show_col_types = FALSE)

CRISPR_FILE <- file.path(DATA_DIR, "CRISPRGeneDependency.csv")
EXPRESSION_FILE <- file.path(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv")
MUTATION_FILE <- file.path(DATA_DIR, "OmicsSomaticMutationsMAF.maf")

has_crispr <- file.exists(CRISPR_FILE)
has_expression <- file.exists(EXPRESSION_FILE)
has_mutation <- file.exists(MUTATION_FILE)

clean_depmap_gene <- function(x) {
  x %>% str_replace(" \\(.*\\)$", "") %>% str_trim()
}

if (has_crispr) {
  CRISPRGeneDependency <- readr::read_csv(
    CRISPR_FILE,
    show_col_types = FALSE,
    name_repair = "unique"
  )
  CRISPR_Dep <- CRISPRGeneDependency
  names(CRISPR_Dep)[1] <- "ModelID"
  names(CRISPR_Dep) <- c("ModelID", clean_depmap_gene(names(CRISPR_Dep)[-1]))
} else {
  CRISPRGeneDependency <- tibble()
  CRISPR_Dep <- tibble()
}

if (has_expression) {
  Expression_TPM <- readr::read_csv(
    EXPRESSION_FILE,
    show_col_types = FALSE,
    name_repair = "unique"
  )

  expression_model_candidates <- names(Expression_TPM)[
    vapply(
      Expression_TPM,
      function(x) any(str_detect(as.character(x), "^ACH-[0-9]+$"), na.rm = TRUE),
      logical(1)
    )
  ]

  if (length(expression_model_candidates) == 0) {
    stop("Could not identify the ACH-* ModelID column in the expression file.")
  }

  expression_model_col <- expression_model_candidates[[1]]
  Expression_TPM <- Expression_TPM %>% rename(DepMapModelID = all_of(expression_model_col))
  names(Expression_TPM) <- clean_depmap_gene(names(Expression_TPM))
} else {
  Expression_TPM <- tibble()
  expression_model_col <- NA_character_
}

if (has_mutation) {
  Mutation_MAF <- read.delim(
    MUTATION_FILE,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>% as_tibble()
} else {
  Mutation_MAF <- tibble()
}

find_first_column <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

mutation_model_col <- if (has_mutation) {
  find_first_column(
    Mutation_MAF,
    c("ModelID", "DepMap_ID", "DepMapID", "Tumor_Sample_Barcode", "SampleID", "ProfileID", "SequencingID", "StrippedCellLineName")
  )
} else NA_character_

mutation_gene_col <- if (has_mutation) {
  find_first_column(Mutation_MAF, c("Hugo_Symbol", "Gene", "GeneSymbol", "gene_symbol"))
} else NA_character_

mutation_class_col <- if (has_mutation) {
  find_first_column(Mutation_MAF, c("Variant_Classification", "Consequence", "VariantType", "Variant_Type"))
} else NA_character_

mutation_protein_col <- if (has_mutation) {
  find_first_column(Mutation_MAF, c("HGVSp_Short", "Protein_Change", "HGVSp", "ProteinChange"))
} else NA_character_

# ============================================================
# 2. CISPLATIN REFERENCE COMBINATIONS
# ============================================================

Clinical_Cisplatin_Combinations_Final <- tibble::tribble(
  ~ClinicalCombination, ~TherapeuticClass, ~DrugNames, ~CompoundIDs,
  "Cisplatin_C225", "EGFR-targeted antibody", c("CISPLATIN", "C225"), c("DPC-001793", "DPC-001418"),
  "Cisplatin_Paclitaxel", "Microtubule inhibitor", c("CISPLATIN", "PACLITAXEL"), c("DPC-001793", "DPC-004880"),
  "Cisplatin_Topotecan", "Topoisomerase I inhibitor", c("CISPLATIN", "TOPOTECAN"), c("DPC-001793", "DPC-006564"),
  "Cisplatin_Fluorouracil", "Antimetabolite", c("CISPLATIN", "FLUOROURACIL"), c("DPC-001793", "DPC-002828"),
  "Cisplatin_Methotrexate", "Antifolate", c("CISPLATIN", "METHOTREXATE"), c("DPC-001793", "DPC-000234"),
  "Cisplatin_Bleomycin", "Cytotoxic chemotherapy", c("CISPLATIN", "BLEOMYCIN"), c("DPC-001793", "DPC-001132"),
  "Cisplatin_Vinorelbine", "Microtubule inhibitor", c("CISPLATIN", "VINORELBINE"), c("DPC-001793", "DPC-002684"),
  "Cisplatin_Vinblastine", "Vinca alkaloid", c("CISPLATIN", "VINBLASTINE"), c("DPC-001793", "DPC-006819"),
  "Cisplatin_Paclitaxel_Veliparib", "PARP inhibitor combination", c("CISPLATIN", "PACLITAXEL", "VELIPARIB"), c("DPC-001793", "DPC-004880", "DPC-000236")
)

# ============================================================
# 3. PREPARE CERVIX DATA
# ============================================================

Cervix_Model <- Model %>%
  filter(OncotreeLineage == "Cervix") %>%
  select(ModelID, CellLineName, StrippedCellLineName, DepmapModelType)

Cervix_GDSC_All <- GDSC2ResponseCurves %>%
  inner_join(Cervix_Model, by = "ModelID") %>%
  distinct(
    CompoundID, CompoundName, SampleID, ModelID, CellLineName,
    EC50, LowerAsymptote, UpperAsymptote, Slope,
    MinimumDose, MaximumDose, DoseUnit,
    .keep_all = TRUE
  )

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
  drug_name_upper <- str_to_upper(drug_name)

  df %>%
    left_join(
      PortalCompounds %>% select(CompoundID, Synonyms, ChEMBLID, PubChemCID, TargetOrMechanism),
      by = "CompoundID"
    ) %>%
    filter(
      case_when(
        drug_name_upper == "CISPLATIN" ~ CompoundID == "DPC-001793",
        drug_name_upper == "C225" ~ CompoundID == "DPC-001418",
        drug_name_upper == "CETUXIMAB" ~ CompoundID == "DPC-001418",
        drug_name_upper == "PACLITAXEL" ~ CompoundID == "DPC-004880",
        drug_name_upper == "TOPOTECAN" ~ CompoundID == "DPC-006564",
        drug_name_upper == "FLUOROURACIL" ~ CompoundID == "DPC-002828",
        drug_name_upper == "METHOTREXATE" ~ CompoundID == "DPC-000234",
        drug_name_upper == "BLEOMYCIN" ~ CompoundID == "DPC-001132",
        drug_name_upper == "VINORELBINE" ~ CompoundID == "DPC-002684",
        drug_name_upper == "VINBLASTINE" ~ CompoundID == "DPC-006819",
        drug_name_upper == "VELIPARIB" ~ CompoundID == "DPC-000236",
        TRUE ~ FALSE
      )
    )
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

# ============================================================
# 5. AVAILABLE THIRD DRUGS AND ANTI-ANGIOGENIC SUBSET
# ============================================================

Available_Third_Drugs <- Cervix_GDSC_All %>%
  distinct(CompoundID, CompoundName) %>%
  left_join(
    PortalCompounds %>%
      select(
        CompoundID,
        Synonyms,
        GeneSymbolOfTargets,
        TargetOrMechanism,
        ChEMBLID,
        PubChemCID
      ),
    by = "CompoundID"
  ) %>%
  arrange(CompoundName)

Third_Drug_Choices <- setNames(
  Available_Third_Drugs$CompoundID,
  paste0(
    Available_Third_Drugs$CompoundName,
    " [",
    Available_Third_Drugs$CompoundID,
    "]"
  )
)

# Terms used to identify VEGF/VEGFR and broader angiogenesis-associated
# compounds in the compound catalogue. The final subset is restricted to
# compounds that have an actual GDSC2 response curve in at least one
# cervical cancer model.
Antiangiogenic_Search_Terms <- c(
  "bevacizumab",
  "avastin",
  "vegf",
  "vegfa",
  "vegfb",
  "vegfc",
  "vegfd",
  "vegfr",
  "vegfr1",
  "vegfr2",
  "vegfr3",
  "flt1",
  "flt4",
  "kdr",
  "vascular endothelial",
  "angiogenesis"
)

make_annotation_text <- function(...) {
  values <- list(...)
  values <- vapply(
    values,
    function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      paste(x, collapse = " ")
    },
    character(1)
  )

  stringr::str_to_lower(paste(values, collapse = " "))
}

identify_matching_terms <- function(search_text, terms) {
  hits <- terms[
    vapply(
      terms,
      function(term) {
        stringr::str_detect(
          search_text,
          stringr::fixed(stringr::str_to_lower(term))
        )
      },
      logical(1)
    )
  ]

  paste(unique(hits), collapse = "; ")
}

Cervical_Antiangiogenic_Compounds <- Available_Third_Drugs %>%
  rowwise() %>%
  mutate(
    AntiangiogenicSearchText = make_annotation_text(
      CompoundName,
      Synonyms,
      GeneSymbolOfTargets,
      TargetOrMechanism
    ),
    MatchedAntiangiogenicAnnotation = identify_matching_terms(
      AntiangiogenicSearchText,
      Antiangiogenic_Search_Terms
    )
  ) %>%
  ungroup() %>%
  filter(MatchedAntiangiogenicAnnotation != "") %>%
  select(-AntiangiogenicSearchText) %>%
  distinct(CompoundID, .keep_all = TRUE) %>%
  arrange(CompoundName)

Antiangiogenic_Drug_Choices <- setNames(
  Cervical_Antiangiogenic_Compounds$CompoundID,
  paste0(
    Cervical_Antiangiogenic_Compounds$CompoundName,
    " [",
    Cervical_Antiangiogenic_Compounds$CompoundID,
    "]"
  )
)

message(
  "Identified ",
  nrow(Cervical_Antiangiogenic_Compounds),
  " anti-angiogenic compounds with cervical GDSC2 response data."
)

# ============================================================
# 6. TARGET-LEVEL VALIDATION MODULES
# ============================================================

extract_targets <- function(target_string) {
  if (length(target_string) == 0 || is.null(target_string) || is.na(target_string[[1]]) || target_string[[1]] == "") {
    return(character(0))
  }

  target_string[[1]] %>%
    str_replace_all(",", ";") %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x == "" || is.na(.x)) %>%
    unique()
}

add_target_annotation <- function(top10_tbl) {
  if (is.null(top10_tbl) || nrow(top10_tbl) == 0) return(tibble())

  top10_tbl %>%
    left_join(
      PortalCompounds %>% select(CompoundID, GeneSymbolOfTargets, TargetOrMechanism),
      by = c("ThirdDrugCompoundID" = "CompoundID")
    ) %>%
    mutate(TargetGene = map(GeneSymbolOfTargets, extract_targets)) %>%
    unnest(TargetGene, keep_empty = TRUE) %>%
    mutate(TargetGene = ifelse(is.na(TargetGene) | TargetGene == "", NA_character_, TargetGene))
}

get_crispr_dependency <- function(model_id, target_gene) {
  if (!has_crispr || nrow(CRISPR_Dep) == 0) {
    return(tibble(CRISPR_Dependency = NA_real_, DependencyClass = "CRISPR file not loaded"))
  }

  if (length(target_gene) == 0 || is.null(target_gene) || is.na(target_gene[[1]]) || target_gene[[1]] == "") {
    return(tibble(CRISPR_Dependency = NA_real_, DependencyClass = "No target annotation"))
  }

  target_gene <- target_gene[[1]]

  if (!target_gene %in% names(CRISPR_Dep)) {
    return(tibble(CRISPR_Dependency = NA_real_, DependencyClass = "Not available"))
  }

  dep_value <- CRISPR_Dep %>% filter(ModelID == model_id) %>% pull(all_of(target_gene))
  dep_value <- suppressWarnings(as.numeric(dep_value))
  dep_value <- if (length(dep_value) == 0) NA_real_ else dep_value[[1]]

  tibble(
    CRISPR_Dependency = dep_value,
    DependencyClass = case_when(
      is.na(CRISPR_Dependency) ~ "Not available",
      CRISPR_Dependency >= 0.8 ~ "Strong dependency",
      CRISPR_Dependency >= 0.5 ~ "Moderate dependency",
      CRISPR_Dependency > 0 ~ "Weak dependency",
      TRUE ~ "Not dependent"
    )
  )
}

get_target_expression <- function(model_id, target_gene) {
  if (!has_expression || nrow(Expression_TPM) == 0) {
    return(tibble(Expression_log2TPM = NA_real_, ExpressionPercentile = NA_real_, ExpressionClass = "Expression file not loaded"))
  }

  if (length(target_gene) == 0 || is.null(target_gene) || is.na(target_gene[[1]]) || target_gene[[1]] == "") {
    return(tibble(Expression_log2TPM = NA_real_, ExpressionPercentile = NA_real_, ExpressionClass = "No target annotation"))
  }

  target_gene <- target_gene[[1]]

  if (!target_gene %in% names(Expression_TPM)) {
    return(tibble(Expression_log2TPM = NA_real_, ExpressionPercentile = NA_real_, ExpressionClass = "Gene not available"))
  }

  model_value <- Expression_TPM %>% filter(DepMapModelID == model_id) %>% pull(all_of(target_gene))
  model_value <- suppressWarnings(as.numeric(model_value))

  if (length(model_value) == 0) {
    return(tibble(Expression_log2TPM = NA_real_, ExpressionPercentile = NA_real_, ExpressionClass = "Model not available"))
  }

  value <- model_value[[1]]

  if (is.na(value) || !is.finite(value)) {
    return(tibble(Expression_log2TPM = NA_real_, ExpressionPercentile = NA_real_, ExpressionClass = "Expression not available"))
  }

  all_values <- suppressWarnings(as.numeric(Expression_TPM[[target_gene]]))
  all_values <- all_values[is.finite(all_values)]
  percentile <- if (length(all_values) > 0) 100 * mean(all_values <= value) else NA_real_

  expr_class <- case_when(
    is.na(percentile) ~ "Not available",
    percentile >= 75 ~ "High expression",
    percentile >= 25 ~ "Moderate expression",
    TRUE ~ "Low expression"
  )

  tibble(Expression_log2TPM = value, ExpressionPercentile = percentile, ExpressionClass = expr_class)
}

get_target_mutation <- function(model_id, target_gene) {
  if (!has_mutation || nrow(Mutation_MAF) == 0) {
    return(tibble(
      MutationStatus = "Mutation file not loaded",
      MutationCount = NA_integer_,
      VariantClassification = NA_character_,
      ProteinChange = NA_character_
    ))
  }

  if (length(target_gene) == 0 || is.null(target_gene) || is.na(target_gene[[1]]) || target_gene[[1]] == "") {
    return(tibble(
      MutationStatus = "No target annotation",
      MutationCount = NA_integer_,
      VariantClassification = NA_character_,
      ProteinChange = NA_character_
    ))
  }

  target_gene <- target_gene[[1]]

  if (is.na(mutation_model_col) || is.na(mutation_gene_col)) {
    return(tibble(
      MutationStatus = "MAF columns not recognised",
      MutationCount = NA_integer_,
      VariantClassification = NA_character_,
      ProteinChange = NA_character_
    ))
  }

  hits <- Mutation_MAF %>%
    filter(
      .data[[mutation_model_col]] == model_id,
      str_to_upper(.data[[mutation_gene_col]]) == str_to_upper(target_gene)
    )

  if (nrow(hits) == 0) {
    return(tibble(
      MutationStatus = "Wild-type/no recorded mutation",
      MutationCount = 0L,
      VariantClassification = NA_character_,
      ProteinChange = NA_character_
    ))
  }

  variant_text <- if (!is.na(mutation_class_col)) {
    hits[[mutation_class_col]] %>% as.character() %>% na.omit() %>% unique() %>% paste(collapse = "; ")
  } else NA_character_

  protein_text <- if (!is.na(mutation_protein_col)) {
    hits[[mutation_protein_col]] %>% as.character() %>% na.omit() %>% unique() %>% paste(collapse = "; ")
  } else NA_character_

  tibble(
    MutationStatus = "Mutated",
    MutationCount = nrow(hits),
    VariantClassification = ifelse(variant_text == "", NA_character_, variant_text),
    ProteinChange = ifelse(protein_text == "", NA_character_, protein_text)
  )
}

build_validation_tables <- function(top10_tbl) {
  target_tbl <- add_target_annotation(top10_tbl)

  if (nrow(target_tbl) == 0) {
    return(list(functional = tibble(), molecular = tibble(), biomarker = tibble()))
  }

  functional_tbl <- target_tbl %>%
    mutate(Functional_Profile = map2(ModelID, TargetGene, get_crispr_dependency)) %>%
    unnest(Functional_Profile)

  molecular_tbl <- functional_tbl %>%
    mutate(Molecular_Profile = map2(ModelID, TargetGene, get_target_expression)) %>%
    unnest(Molecular_Profile)

  biomarker_tbl <- molecular_tbl %>%
    mutate(Biomarker_Profile = map2(ModelID, TargetGene, get_target_mutation)) %>%
    unnest(Biomarker_Profile)

  list(functional = functional_tbl, molecular = molecular_tbl, biomarker = biomarker_tbl)
}

# ============================================================
# 7. SHINY UI
# ============================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Cisplatin Clinical Reference Combination Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput(
        "combo",
        "Select clinical cisplatin reference combination:",
        choices = Clinical_Cisplatin_Combinations_Final$ClinicalCombination,
        selected = "Cisplatin_Paclitaxel"
      ),

      hr(),

      radioButtons(
        "third_drug_mode",
        "Third-drug screening mode:",
        choices = c(
          "Screen all cervical GDSC2 compounds" = "all",
          "Screen all anti-angiogenic cervical compounds" = "antiangiogenic",
          "Select third drugs manually" = "selected"
        ),
        selected = "antiangiogenic"
      ),

      conditionalPanel(
        condition = "input.third_drug_mode == 'antiangiogenic'",
        wellPanel(
          tags$b("Anti-angiogenic screening set"),
          tags$p(
            paste0(
              nrow(Cervical_Antiangiogenic_Compounds),
              " compounds with VEGF-, VEGFR-, KDR-, FLT1-, FLT4- or angiogenesis-related annotations and cervical GDSC2 response data will be selected automatically."
            )
          ),
          tags$p(
            "The subset is generated dynamically from PortalCompounds.csv and Cervix_GDSC_All, so only compounds tested in at least one cervical cancer model are included."
          )
        )
      ),

      conditionalPanel(
        condition = "input.third_drug_mode == 'selected'",
        selectizeInput(
          "selected_third_drugs",
          "Search and select third drugs:",
          choices = Third_Drug_Choices,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Type drug name or CompoundID...",
            maxOptions = 1000,
            plugins = list("remove_button")
          )
        ),
        helpText(
          "Each selected drug is evaluated independently as monotherapy and as an addition to the selected clinical reference backbone."
        )
      ),

      hr(),

      selectInput(
        "synergy_method",
        "Synergy model:",
        choices = c(
          "Bliss independence" = "bliss",
          "Highest Single Agent (HSA)" = "hsa",
          "Loewe additivity" = "loewe",
          "ZIP approximation" = "zip"
        ),
        selected = "bliss"
      ),

      helpText(
        paste(
          "The selected model defines the expected non-interacting response.",
          "Synergy score = observed response - model-specific expected response.",
          "The ZIP option is a transparent curve-based approximation because",
          "the app reconstructs response surfaces from fitted single-agent curves."
        )
      ),

      hr(),

      checkboxGroupInput(
        "doses",
        "Dose levels:",
        choices = c(0.001, 0.01, 0.1, 1),
        selected = c(0.001, 0.01, 0.1, 1)
      ),

      numericInput("noise", "Noise SD:", value = 5, min = 0, max = 20, step = 1),
      actionButton("run", "Run analysis", class = "btn-primary"),

      br(), br(),
      downloadButton(
        "download_antiangiogenic_set",
        "Download anti-angiogenic compound set"
      ),
      br(), br(),
      downloadButton("download_reference", "Download reference summary"),
      br(), br(),
      downloadButton("download_monotherapy", "Download monotherapy summary"),
      br(), br(),
      downloadButton("download_monotherapy_dose", "Download monotherapy dose-response"),
      br(), br(),
      downloadButton("download_ranking", "Download third-drug ranking"),
      br(), br(),
      downloadButton("download_top10", "Download top 10"),
      br(), br(),
      downloadButton("download_functional", "Download functional CRISPR"),
      br(), br(),
      downloadButton("download_molecular", "Download molecular expression"),
      br(), br(),
      downloadButton("download_biomarker", "Download optional biomarker")
    ),

    mainPanel(
      verbatimTextOutput("run_status"),

      tabsetPanel(
        tabPanel("Eligible models", DTOutput("eligible_table")),
        tabPanel("Reference summary", DTOutput("reference_table")),
        tabPanel("Selected third drugs", DTOutput("third_drug_selection_table")),

        tabPanel(
          "Monotherapy response",
          h4("How strongly does each third drug inhibit the matched cervical cancer model by itself?"),
          DTOutput("monotherapy_table"),
          plotOutput("monotherapy_plot", height = "850px")
        ),

        tabPanel(
          "Monotherapy dose-response",
          h4("Predicted single-drug inhibition across the selected dose levels"),
          DTOutput("monotherapy_dose_table"),
          plotOutput("monotherapy_dose_plot", height = "850px")
        ),

        tabPanel("Third-drug ranking", DTOutput("ranking_table")),
        tabPanel("Top candidates per model", DTOutput("top10_table")),
        tabPanel("Synergy plot", plotOutput("top10_plot", height = "850px")),

        tabPanel(
          "Mono vs combination",
          h4("Comparison of intrinsic monotherapy response with observed combination response"),
          plotOutput("mono_combination_plot", height = "850px")
        ),

        tabPanel(
          "Functional: CRISPR dependency",
          h4("Does the model require the drug target for survival?"),
          DTOutput("functional_table"),
          plotOutput("functional_plot", height = "750px")
        ),

        tabPanel(
          "Molecular: target expression",
          h4("Is the target actually expressed in that model?"),
          DTOutput("molecular_table"),
          plotOutput("molecular_plot", height = "750px")
        ),

        tabPanel(
          "Optional biomarker: mutations",
          h4("Is dependency or synergy improvement associated with a genomic alteration?"),
          DTOutput("biomarker_table"),
          plotOutput("biomarker_plot", height = "650px")
        )
      )
    )
  )
)

# ============================================================
# 8. SHINY SERVER
# ============================================================

server <- function(input, output, session) {
  selected_combo <- reactive({
    Clinical_Cisplatin_Combinations_Final %>%
      filter(ClinicalCombination == input$combo)
  })

  selected_third_drugs <- reactive({
    if (input$third_drug_mode == "all") {
      return(
        Available_Third_Drugs %>%
          mutate(
            SelectionStatus =
              "All available cervical GDSC2 compounds"
          )
      )
    }

    if (input$third_drug_mode == "antiangiogenic") {
      return(
        Cervical_Antiangiogenic_Compounds %>%
          mutate(
            SelectionStatus =
              "Automatically selected anti-angiogenic cervical compound"
          )
      )
    }

    if (input$third_drug_mode == "selected") {
      req(input$selected_third_drugs)

      return(
        Available_Third_Drugs %>%
          filter(CompoundID %in% input$selected_third_drugs) %>%
          mutate(SelectionStatus = "User selected")
      )
    }

    tibble()
  })

  results <- eventReactive(input$run, {
    validate(need(length(input$doses) > 0, "Please select at least one dose level."))

    combo_row <- selected_combo()
    selected_drugs <- selected_third_drugs()

    if (input$third_drug_mode %in% c("antiangiogenic", "selected")) {
      empty_message <- if (input$third_drug_mode == "antiangiogenic") {
        paste(
          "No anti-angiogenic compounds with cervical GDSC2 response",
          "data were identified. Check the compound annotations and files."
        )
      } else {
        "Please select at least one third drug."
      }

      validate(need(nrow(selected_drugs) > 0, empty_message))
      selected_partner_ids <- unique(selected_drugs$CompoundID)
    } else {
      selected_partner_ids <- NULL
    }

    withProgress(message = "Running pharmacology and validation layers...", value = 0, {
      incProgress(0.20, detail = paste("Running reference-regimen", synergy_method_label(input$synergy_method), "analysis"))

      out <- run_reference_pipeline(
        combo_name = combo_row$ClinicalCombination[[1]],
        therapeutic_class = combo_row$TherapeuticClass[[1]],
        drug_names = combo_row$DrugNames[[1]],
        compound_ids = combo_row$CompoundIDs[[1]],
        doses = as.numeric(input$doses),
        noise_sd = input$noise,
        synergy_method = input$synergy_method,
        selected_partner_ids = selected_partner_ids
      )

      incProgress(0.20, detail = "Evaluating monotherapy response")
      out$third_drug_selection <- selected_drugs

      incProgress(0.20, detail = "Building target-level validation tables")
      validation <- build_validation_tables(out$top10)

      incProgress(0.15, detail = "Adding CRISPR functional evidence")
      out$functional_tbl <- validation$functional

      incProgress(0.15, detail = "Adding target-expression evidence")
      out$molecular_tbl <- validation$molecular

      incProgress(0.10, detail = "Adding mutation biomarker context")
      out$biomarker_tbl <- validation$biomarker

      out
    })
  })

  output$run_status <- renderText({
    if (is.null(results())) return("Click 'Run analysis' to start.")

    out <- results()
    screening_mode <- switch(
      input$third_drug_mode,
      all = "All cervical GDSC2 compounds",
      antiangiogenic = paste0(
        "All anti-angiogenic cervical compounds (n = ",
        nrow(Cervical_Antiangiogenic_Compounds),
        ")"
      ),
      selected = "User-selected third drugs"
    )

    paste(
      paste0("Screening mode: ", screening_mode),
      paste0("Synergy model: ", synergy_method_label(input$synergy_method)),
      paste0("Selected third drugs: ", nrow(out$third_drug_selection)),
      paste0("Eligible models: ", nrow(out$eligible_models)),
      paste0("Reference summary rows: ", nrow(out$reference_summary)),
      paste0("Monotherapy summary rows: ", nrow(out$monotherapy_summary)),
      paste0("Monotherapy dose-response rows: ", nrow(out$monotherapy_dose_response)),
      paste0("Third-drug ranking rows: ", nrow(out$third_ranking)),
      paste0("Top candidate rows: ", nrow(out$top10)),
      paste0("Functional CRISPR rows: ", nrow(out$functional_tbl)),
      paste0("Molecular expression rows: ", nrow(out$molecular_tbl)),
      paste0("Optional biomarker rows: ", nrow(out$biomarker_tbl)),
      paste0("CRISPR loaded: ", has_crispr),
      paste0("Expression loaded: ", has_expression),
      paste0("Mutation MAF loaded: ", has_mutation),
      sep = "\n"
    )
  })

  safe_dt <- function(dat, msg, page_length = 25) {
    if (is.null(dat) || nrow(dat) == 0) dat <- tibble(Message = msg)
    datatable(dat, filter = "top", options = list(pageLength = page_length, scrollX = TRUE))
  }

  output$eligible_table <- renderDT({ req(results()); safe_dt(results()$eligible_models, "No eligible models.", 10) })
  output$reference_table <- renderDT({ req(results()); safe_dt(results()$reference_summary, "No reference summary available.", 10) })
  output$third_drug_selection_table <- renderDT({ req(results()); safe_dt(results()$third_drug_selection, "No third drugs selected.", 25) })
  output$monotherapy_table <- renderDT({ req(results()); safe_dt(results()$monotherapy_summary, "No monotherapy response data available.", 25) })
  output$monotherapy_dose_table <- renderDT({ req(results()); safe_dt(results()$monotherapy_dose_response, "No monotherapy dose-response data available.", 25) })
  output$ranking_table <- renderDT({ req(results()); safe_dt(results()$third_ranking, "No third-drug ranking available.", 25) })
  output$top10_table <- renderDT({ req(results()); safe_dt(results()$top10, "No candidate drugs available.", 25) })
  output$functional_table <- renderDT({ req(results()); safe_dt(results()$functional_tbl, "No functional CRISPR dependency evidence available.", 25) })
  output$molecular_table <- renderDT({ req(results()); safe_dt(results()$molecular_tbl, "No molecular target-expression evidence available.", 25) })
  output$biomarker_table <- renderDT({ req(results()); safe_dt(results()$biomarker_tbl, "No optional mutation biomarker evidence available.", 25) })

  output$monotherapy_plot <- renderPlot({
    req(results())
    plot_data <- results()$monotherapy_summary
    validate(need(nrow(plot_data) > 0, "No eligible monotherapy results are available for plotting."))

    plot_data <- plot_data %>%
      group_by(CellLineName) %>%
      slice_max(MeanMonotherapyResponse, n = 10, with_ties = FALSE) %>%
      ungroup()

    ggplot(plot_data, aes(x = reorder(ThirdDrug, MeanMonotherapyResponse), y = MeanMonotherapyResponse, fill = MonotherapyClass)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Top monotherapy responses in matched cervical models",
        subtitle = "Predicted inhibition is derived from the fitted GDSC2 response curve",
        x = "Third drug",
        y = "Mean predicted monotherapy inhibition (%)",
        fill = "Monotherapy class"
      )
  })

  output$monotherapy_dose_plot <- renderPlot({
    req(results())
    plot_data <- results()$monotherapy_dose_response
    validate(need(nrow(plot_data) > 0, "No monotherapy dose-response data are available for plotting."))

    top_ids <- results()$monotherapy_summary %>%
      group_by(CellLineName) %>%
      slice_max(MeanMonotherapyResponse, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      select(ModelID, ThirdDrugCompoundID)

    plot_data <- plot_data %>%
      semi_join(top_ids, by = c("ModelID", "ThirdDrugCompoundID"))

    ggplot(plot_data, aes(x = Dose, y = PredictedMonotherapyInhibition, group = ThirdDrug, linetype = ThirdDrug)) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.6) +
      scale_x_log10() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = "Third-drug monotherapy dose-response profiles",
        subtitle = "Only the top ten monotherapy candidates per model are shown",
        x = "Dose",
        y = "Predicted monotherapy inhibition (%)",
        linetype = "Third drug"
      )
  })

  output$top10_plot <- renderPlot({
    req(results())
    plot_data <- results()$top10
    validate(need(nrow(plot_data) > 0, "No eligible third-drug results are available for plotting."))

    reference_lines <- results()$reference_summary %>%
      distinct(ModelID, CellLineName, ReferenceMeanSynergy)

    ggplot(plot_data, aes(x = reorder(ThirdDrug, MeanSynergy), y = MeanSynergy, fill = Interaction)) +
      geom_col() +
      geom_hline(
        data = reference_lines,
        aes(yintercept = ReferenceMeanSynergy),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0(
          "Cervical models: ",
          str_replace_all(input$combo, "_", " "),
          " + top third drugs (",
          synergy_method_label(input$synergy_method),
          ")"
        ),
        subtitle = paste(
          "Dashed line = model-specific reference mean;",
          "interaction model =", synergy_method_label(input$synergy_method)
        ),
        x = "Third drug",
        y = "Mean synergy score",
        fill = "Interaction"
      )
  })

  output$mono_combination_plot <- renderPlot({
    req(results())
    plot_data <- results()$top10
    validate(need(nrow(plot_data) > 0, "No monotherapy-versus-combination results are available."))

    ggplot(plot_data, aes(x = MeanMonotherapyResponse, y = MeanObserved, label = ThirdDrug, shape = Interaction)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_point(size = 2.4) +
      geom_text(check_overlap = TRUE, nudge_y = 2, size = 3) +
      facet_wrap(~ CellLineName, scales = "free") +
      theme_bw(base_size = 11) +
      labs(
        title = "Intrinsic monotherapy response versus observed combination response",
        subtitle = "Points above the diagonal show greater observed combination inhibition than monotherapy inhibition",
        x = "Mean predicted monotherapy inhibition (%)",
        y = "Mean observed combination inhibition (%)",
        shape = "Interaction"
      )
  })

  output$functional_plot <- renderPlot({
    req(results())
    plot_data <- results()$functional_tbl %>%
      filter(!is.na(CRISPR_Dependency), DependencyClass %in% c("Strong dependency", "Moderate dependency", "Weak dependency"))

    validate(need(nrow(plot_data) > 0, "No CRISPR dependency data available for plotted targets."))

    ggplot(plot_data, aes(x = reorder(TargetGene, CRISPR_Dependency), y = CRISPR_Dependency, fill = DependencyClass)) +
      geom_col() +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 0.8, linetype = "dotted") +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Functional evidence: CRISPR dependency of third-drug targets",
        subtitle = "Dashed = moderate dependency cutoff; dotted = strong dependency cutoff",
        x = "Target gene",
        y = "CRISPR dependency score",
        fill = "Dependency class"
      )
  })

  output$molecular_plot <- renderPlot({
    req(results())
    plot_data <- results()$molecular_tbl %>% filter(!is.na(Expression_log2TPM))
    validate(need(nrow(plot_data) > 0, "No expression data available for plotted targets."))

    ggplot(plot_data, aes(x = reorder(TargetGene, Expression_log2TPM), y = Expression_log2TPM, fill = ExpressionClass)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Molecular evidence: target expression in matched cervical models",
        subtitle = "Expression values are log2(TPM + 1); expression classes use target-specific percentiles",
        x = "Target gene",
        y = "log2(TPM + 1)",
        fill = "Expression class"
      )
  })

  output$biomarker_plot <- renderPlot({
    req(results())
    plot_data <- results()$biomarker_tbl %>% count(CellLineName, MutationStatus, name = "n")
    validate(need(nrow(plot_data) > 0, "No mutation biomarker data available."))

    ggplot(plot_data, aes(x = CellLineName, y = n, fill = MutationStatus)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 11) +
      labs(
        title = "Optional biomarker evidence: target mutation context",
        x = "Cervical cancer model",
        y = "Drug-target associations",
        fill = "Mutation status"
      )
  })

  output$download_antiangiogenic_set <- downloadHandler(
    filename = function() {
      "Cervical_GDSC2_Antiangiogenic_Compound_Set.csv"
    },
    content = function(file) {
      write.csv(
        Cervical_Antiangiogenic_Compounds,
        file,
        row.names = FALSE
      )
    }
  )

  output$download_reference <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_Reference_Summary.csv"
    ),
    content = function(file) write.csv(results()$reference_summary, file, row.names = FALSE)
  )

  output$download_monotherapy <- downloadHandler(
    filename = function() paste0(input$combo, "_Monotherapy_Response_Summary.csv"),
    content = function(file) write.csv(results()$monotherapy_summary, file, row.names = FALSE)
  )

  output$download_monotherapy_dose <- downloadHandler(
    filename = function() paste0(input$combo, "_Monotherapy_Dose_Response.csv"),
    content = function(file) write.csv(results()$monotherapy_dose_response, file, row.names = FALSE)
  )

  output$download_ranking <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_ThirdDrug_Ranking.csv"
    ),
    content = function(file) write.csv(results()$third_ranking, file, row.names = FALSE)
  )

  output$download_top10 <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_Top_Candidates_Per_Model.csv"
    ),
    content = function(file) write.csv(results()$top10, file, row.names = FALSE)
  )

  output$download_functional <- downloadHandler(
    filename = function() paste0(input$combo, "_Functional_CRISPR_Evidence.csv"),
    content = function(file) write.csv(results()$functional_tbl, file, row.names = FALSE)
  )

  output$download_molecular <- downloadHandler(
    filename = function() paste0(input$combo, "_Molecular_Expression_Evidence.csv"),
    content = function(file) write.csv(results()$molecular_tbl, file, row.names = FALSE)
  )

  output$download_biomarker <- downloadHandler(
    filename = function() paste0(input$combo, "_Optional_Mutation_Biomarker_Evidence.csv"),
    content = function(file) write.csv(results()$biomarker_tbl, file, row.names = FALSE)
  )
}

# ============================================================
# 9. RUN APP
# ============================================================

shinyApp(ui = ui, server = server)
