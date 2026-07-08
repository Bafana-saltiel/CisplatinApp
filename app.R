# ============================================================
# Shiny App: Cisplatin Clinical Reference Combination Explorer
# Updated validation-layer version
#
# Validation panels:
# 1. Functional evidence: CRISPRGeneDependency.csv
# 2. Molecular evidence: OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv
# 3. Optional biomarker evidence: OmicsSomaticMutationsMAF.maf
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
## ============================================================
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

PortalCompounds <- readr::read_csv(
  required_files[["PortalCompounds"]],
  show_col_types = FALSE
)

Model <- readr::read_csv(
  required_files[["Model"]],
  show_col_types = FALSE
)

GDSC2ResponseCurves <- readr::read_csv(
  required_files[["GDSC2ResponseCurves"]],
  show_col_types = FALSE
)

CRISPR_FILE <- file.path(DATA_DIR, "CRISPRGeneDependency.csv")

EXPRESSION_FILE <- file.path(
  DATA_DIR,
  "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"
)

MUTATION_FILE <- file.path(
  DATA_DIR,
  "OmicsSomaticMutationsMAF.maf"
)

has_crispr <- file.exists(CRISPR_FILE)
has_expression <- file.exists(EXPRESSION_FILE)
has_mutation <- file.exists(MUTATION_FILE)

clean_depmap_gene <- function(x) {
  x %>%
    stringr::str_replace(" \\(.*\\)$", "") %>%
    stringr::str_trim()
}

# -----------------------------
# CRISPR dependency
# -----------------------------

if (has_crispr) {
  
  CRISPRGeneDependency <- readr::read_csv(
    CRISPR_FILE,
    show_col_types = FALSE,
    name_repair = "unique"
  )
  
  CRISPR_Dep <- CRISPRGeneDependency
  
  names(CRISPR_Dep)[1] <- "ModelID"
  
  names(CRISPR_Dep) <- c(
    "ModelID",
    clean_depmap_gene(names(CRISPR_Dep)[-1])
  )
  
} else {
  
  CRISPRGeneDependency <- tibble::tibble()
  CRISPR_Dep <- tibble::tibble()
}

# -----------------------------
# Expression data
# -----------------------------

if (has_expression) {
  
  Expression_TPM <- readr::read_csv(
    EXPRESSION_FILE,
    show_col_types = FALSE,
    name_repair = "unique"
  )
  
  expression_model_candidates <- names(Expression_TPM)[
    vapply(
      Expression_TPM,
      function(x) {
        any(
          stringr::str_detect(
            as.character(x),
            "^ACH-[0-9]+$"
          ),
          na.rm = TRUE
        )
      },
      logical(1)
    )
  ]
  
  if (length(expression_model_candidates) == 0) {
    stop(
      "Could not identify the ACH-* ModelID column in the expression file."
    )
  }
  
  expression_model_col <- expression_model_candidates[[1]]
  
  Expression_TPM <- Expression_TPM %>%
    dplyr::rename(
      DepMapModelID = dplyr::all_of(expression_model_col)
    )
  
  names(Expression_TPM) <- clean_depmap_gene(names(Expression_TPM))
  
} else {
  
  Expression_TPM <- tibble::tibble()
  expression_model_col <- NA_character_
}

# -----------------------------
# Mutation data
# -----------------------------

if (has_mutation) {
  
  Mutation_MAF <- read.delim(
    MUTATION_FILE,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    tibble::as_tibble()
  
} else {
  
  Mutation_MAF <- tibble::tibble()
}

find_first_column <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    NA_character_
  } else {
    hit[1]
  }
}

mutation_model_col <- if (has_mutation) {
  find_first_column(
    Mutation_MAF,
    c(
      "ModelID",
      "DepMap_ID",
      "DepMapID",
      "Tumor_Sample_Barcode",
      "SampleID",
      "ProfileID",
      "SequencingID",
      "StrippedCellLineName"
    )
  )
} else {
  NA_character_
}

mutation_gene_col <- if (has_mutation) {
  find_first_column(
    Mutation_MAF,
    c(
      "Hugo_Symbol",
      "Gene",
      "GeneSymbol",
      "gene_symbol"
    )
  )
} else {
  NA_character_
}

mutation_class_col <- if (has_mutation) {
  find_first_column(
    Mutation_MAF,
    c(
      "Variant_Classification",
      "Consequence",
      "VariantType",
      "Variant_Type"
    )
  )
} else {
  NA_character_
}

mutation_protein_col <- if (has_mutation) {
  find_first_column(
    Mutation_MAF,
    c(
      "HGVSp_Short",
      "Protein_Change",
      "HGVSp",
      "ProteinChange"
    )
  )
} else {
  NA_character_
}

# ============================================================
# 2. CISPLATIN REFERENCE COMBINATIONS
# ============================================================

Clinical_Cisplatin_Combinations_Final <- tibble::tribble(
  ~ClinicalCombination, ~TherapeuticClass, ~DrugNames, ~CompoundIDs,

  "Cisplatin_C225",
  "EGFR-targeted antibody",
  c("CISPLATIN", "C225"),
  c("DPC-001793", "DPC-001418"),

  "Cisplatin_Paclitaxel",
  "Microtubule inhibitor",
  c("CISPLATIN", "PACLITAXEL"),
  c("DPC-001793", "DPC-004880"),

  "Cisplatin_Topotecan",
  "Topoisomerase I inhibitor",
  c("CISPLATIN", "TOPOTECAN"),
  c("DPC-001793", "DPC-006564"),

  "Cisplatin_Fluorouracil",
  "Antimetabolite",
  c("CISPLATIN", "FLUOROURACIL"),
  c("DPC-001793", "DPC-002828"),

  "Cisplatin_Methotrexate",
  "Antifolate",
  c("CISPLATIN", "METHOTREXATE"),
  c("DPC-001793", "DPC-000234"),

  "Cisplatin_Bleomycin",
  "Cytotoxic chemotherapy",
  c("CISPLATIN", "BLEOMYCIN"),
  c("DPC-001793", "DPC-001132"),

  "Cisplatin_Vinorelbine",
  "Microtubule inhibitor",
  c("CISPLATIN", "VINORELBINE"),
  c("DPC-001793", "DPC-002684"),

  "Cisplatin_Vinblastine",
  "Vinca alkaloid",
  c("CISPLATIN", "VINBLASTINE"),
  c("DPC-001793", "DPC-006819"),

  "Cisplatin_Paclitaxel_Veliparib",
  "PARP inhibitor combination",
  c("CISPLATIN", "PACLITAXEL", "VELIPARIB"),
  c("DPC-001793", "DPC-004880", "DPC-000236")
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
  E <- E / 100
  100 * (1 - prod(1 - E))
}

find_drug_rows <- function(drug_name, df = Cervix_GDSC_All) {

  drug_name_upper <- str_to_upper(drug_name)

  df %>%
    left_join(
      PortalCompounds %>%
        select(
          CompoundID,
          Synonyms,
          ChEMBLID,
          PubChemCID,
          TargetOrMechanism
        ),
      by = "CompoundID"
    ) %>%
    filter(
      case_when(
        drug_name_upper == "CISPLATIN"    ~ CompoundID == "DPC-001793",
        drug_name_upper == "C225"         ~ CompoundID == "DPC-001418",
        drug_name_upper == "CETUXIMAB"    ~ CompoundID == "DPC-001418",
        drug_name_upper == "PACLITAXEL"   ~ CompoundID == "DPC-004880",
        drug_name_upper == "TOPOTECAN"    ~ CompoundID == "DPC-006564",
        drug_name_upper == "FLUOROURACIL" ~ CompoundID == "DPC-002828",
        drug_name_upper == "METHOTREXATE" ~ CompoundID == "DPC-000234",
        drug_name_upper == "BLEOMYCIN"    ~ CompoundID == "DPC-001132",
        drug_name_upper == "VINORELBINE"  ~ CompoundID == "DPC-002684",
        drug_name_upper == "VINBLASTINE"  ~ CompoundID == "DPC-006819",
        drug_name_upper == "VELIPARIB"    ~ CompoundID == "DPC-000236",
        TRUE ~ FALSE
      )
    )
}

simulate_reference_combo <- function(model_id, combo_name, therapeutic_class,
                                     drug_names, doses, noise_sd) {

  df_model <- Cervix_GDSC_All %>%
    filter(ModelID == model_id)

  matched_drugs <- map_dfr(drug_names, function(drug) {
    rows <- find_drug_rows(drug, df_model)
    if (nrow(rows) == 0) return(tibble())
    rows %>% slice(1) %>% mutate(QueryDrug = drug)
  })

  if (nrow(matched_drugs) != length(drug_names)) return(tibble())

  dose_grid <- expand.grid(
    replicate(length(drug_names), doses, simplify = FALSE)
  )

  colnames(dose_grid) <- paste0(matched_drugs$QueryDrug, "_Dose")

  inhibition_matrix <- map_dfc(seq_len(nrow(matched_drugs)), function(i) {
    inhib <- predict_inhibition(
      dose  = dose_grid[[i]],
      ec50  = matched_drugs$EC50[i],
      lower = matched_drugs$LowerAsymptote[i],
      upper = matched_drugs$UpperAsymptote[i],
      slope = matched_drugs$Slope[i]
    )
    tibble(!!paste0(matched_drugs$QueryDrug[i], "_Response") := inhib)
  })

  E_bliss <- apply(inhibition_matrix, 1, bliss_multi)

  set.seed(make_seed(model_id, combo_name, "REFERENCE"))

  E_observed <- pmax(
    pmin(E_bliss + rnorm(length(E_bliss), mean = 0, sd = noise_sd), 100),
    0
  )

  bind_cols(
    tibble(
      ModelID = model_id,
      CellLineName = matched_drugs$CellLineName[1],
      DepmapModelType = matched_drugs$DepmapModelType[1],
      ClinicalCombination = combo_name,
      TherapeuticClass = therapeutic_class,
      RequiredDrugs = paste(drug_names, collapse = " + "),
      MatchedCompoundNames = paste(matched_drugs$CompoundName, collapse = " + "),
      MatchedCompoundIDs = paste(matched_drugs$CompoundID, collapse = " + ")
    ),
    dose_grid,
    inhibition_matrix,
    tibble(
      BlissExpected = E_bliss,
      Observed = E_observed,
      BlissScore = E_observed - E_bliss
    )
  )
}

simulate_reference_plus_partner <- function(model_id, combo_name, therapeutic_class,
                                            drug_names, partner, doses, noise_sd) {

  df_model <- Cervix_GDSC_All %>%
    filter(ModelID == model_id)

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

  dose_grid <- expand.grid(
    replicate(nrow(matched_drugs), doses, simplify = FALSE)
  )

  colnames(dose_grid) <- paste0(matched_drugs$QueryDrug, "_Dose")

  inhibition_matrix <- map_dfc(seq_len(nrow(matched_drugs)), function(i) {
    inhib <- predict_inhibition(
      dose  = dose_grid[[i]],
      ec50  = matched_drugs$EC50[i],
      lower = matched_drugs$LowerAsymptote[i],
      upper = matched_drugs$UpperAsymptote[i],
      slope = matched_drugs$Slope[i]
    )
    tibble(!!paste0(matched_drugs$QueryDrug[i], "_Response") := inhib)
  })

  E_bliss <- apply(inhibition_matrix, 1, bliss_multi)

  set.seed(make_seed(model_id, combo_name, partner$CompoundID[1]))

  E_observed <- pmax(
    pmin(E_bliss + rnorm(length(E_bliss), mean = 0, sd = noise_sd), 100),
    0
  )

  bind_cols(
    tibble(
      ModelID = model_id,
      CellLineName = matched_reference$CellLineName[1],
      DepmapModelType = matched_reference$DepmapModelType[1],
      ReferenceCombination = combo_name,
      TherapeuticClass = therapeutic_class,
      ThirdDrug = partner$CompoundName[1],
      ThirdDrugCompoundID = partner$CompoundID[1],
      ReferenceCompoundNames = paste(matched_reference$CompoundName, collapse = " + "),
      ReferenceCompoundIDs = paste(matched_reference$CompoundID, collapse = " + ")
    ),
    dose_grid,
    inhibition_matrix,
    tibble(
      BlissExpected = E_bliss,
      Observed = E_observed,
      BlissScore = E_observed - E_bliss
    )
  )
}

run_reference_pipeline <- function(combo_name, drug_names, compound_ids,
                                   therapeutic_class, doses, noise_sd,
                                   selected_partner_ids = NULL) {

  eligible_models <- Cervix_GDSC_All %>%
    filter(CompoundID %in% compound_ids) %>%
    distinct(ModelID, CellLineName, DepmapModelType, CompoundID) %>%
    count(ModelID, CellLineName, DepmapModelType, name = "n_reference_drugs") %>%
    filter(n_reference_drugs == length(unique(compound_ids))) %>%
    select(ModelID, CellLineName, DepmapModelType)

  if (nrow(eligible_models) == 0) {
    return(list(
      eligible_models = tibble(),
      reference_bliss = tibble(),
      reference_summary = tibble(),
      third_ranking = tibble(),
      top10 = tibble()
    ))
  }

  reference_bliss <- map_dfr(
    eligible_models$ModelID,
    ~ simulate_reference_combo(
      model_id = .x,
      combo_name = combo_name,
      therapeutic_class = therapeutic_class,
      drug_names = drug_names,
      doses = doses,
      noise_sd = noise_sd
    )
  )

  reference_summary <- reference_bliss %>%
    group_by(
      ModelID, CellLineName, DepmapModelType,
      ClinicalCombination, TherapeuticClass,
      RequiredDrugs, MatchedCompoundNames, MatchedCompoundIDs
    ) %>%
    summarise(
      ReferenceMeanBliss = mean(BlissScore, na.rm = TRUE),
      ReferenceMaxBliss = max(BlissScore, na.rm = TRUE),
      ReferenceMinBliss = min(BlissScore, na.rm = TRUE),
      ReferenceSDBliss = sd(BlissScore, na.rm = TRUE),
      ReferenceMeanExpected = mean(BlissExpected, na.rm = TRUE),
      ReferenceMeanObserved = mean(Observed, na.rm = TRUE),
      .groups = "drop"
    )

  run_model <- function(model_id) {

    df_model <- Cervix_GDSC_All %>%
      filter(ModelID == model_id)

    partners <- df_model %>%
      filter(!CompoundID %in% compound_ids) %>%
      distinct(
        CompoundID, CompoundName, ModelID, CellLineName, DepmapModelType,
        EC50, LowerAsymptote, UpperAsymptote, Slope,
        MinimumDose, MaximumDose, DoseUnit,
        .keep_all = TRUE
      ) %>%
      filter(
        !is.na(EC50),
        !is.na(LowerAsymptote),
        !is.na(UpperAsymptote),
        !is.na(Slope)
      )

    if (!is.null(selected_partner_ids)) {
      partners <- partners %>%
        filter(CompoundID %in% selected_partner_ids)
    }

    if (nrow(partners) == 0) return(tibble())

    map_dfr(
      seq_len(nrow(partners)),
      ~ simulate_reference_plus_partner(
        model_id = model_id,
        combo_name = combo_name,
        therapeutic_class = therapeutic_class,
        drug_names = drug_names,
        partner = partners[.x, ],
        doses = doses,
        noise_sd = noise_sd
      )
    )
  }

  third_bliss <- map_dfr(eligible_models$ModelID, run_model)

  if (nrow(third_bliss) == 0) {
    return(list(
      eligible_models = eligible_models,
      reference_bliss = reference_bliss,
      reference_summary = reference_summary,
      third_ranking = tibble(),
      top10 = tibble()
    ))
  }

  third_ranking <- third_bliss %>%
    group_by(
      ModelID, CellLineName, DepmapModelType,
      ReferenceCombination, TherapeuticClass,
      ThirdDrug, ThirdDrugCompoundID,
      ReferenceCompoundNames, ReferenceCompoundIDs
    ) %>%
    summarise(
      n_dose_combinations = n(),
      MeanBliss = mean(BlissScore, na.rm = TRUE),
      MaxBliss = max(BlissScore, na.rm = TRUE),
      MinBliss = min(BlissScore, na.rm = TRUE),
      SDBliss = sd(BlissScore, na.rm = TRUE),
      MeanExpected = mean(BlissExpected, na.rm = TRUE),
      MeanObserved = mean(Observed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      reference_summary %>% select(ModelID, ReferenceMeanBliss),
      by = "ModelID"
    ) %>%
    mutate(
      Delta_vs_Reference = MeanBliss - ReferenceMeanBliss,
      Interaction = case_when(
        MeanBliss > ReferenceMeanBliss ~ "Stronger_than_reference",
        MeanBliss < 0 ~ "Antagonistic",
        TRUE ~ "Weaker_or_Additive"
      )
    ) %>%
    arrange(CellLineName, desc(MeanBliss))

  top10 <- third_ranking %>%
    group_by(CellLineName) %>%
    slice_max(MeanBliss, n = 10, with_ties = FALSE) %>%
    ungroup()

  list(
    eligible_models = eligible_models,
    reference_bliss = reference_bliss,
    reference_summary = reference_summary,
    third_ranking = third_ranking,
    top10 = top10
  )
}

# ============================================================
# 5. AVAILABLE THIRD DRUGS FOR MULTIPLE SELECTION
# ============================================================

Available_Third_Drugs <- Cervix_GDSC_All %>%
  distinct(CompoundID, CompoundName) %>%
  left_join(
    PortalCompounds %>%
      select(
        CompoundID,
        Synonyms,
        GeneSymbolOfTargets,
        TargetOrMechanism
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

# ============================================================
# 6. TARGET-LEVEL VALIDATION MODULES
# ============================================================

extract_targets <- function(target_string) {
  if (is.na(target_string) || target_string == "") return(character(0))

  target_string %>%
    str_replace_all(",", ";") %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x == "" | is.na(.x)) %>%
    unique()
}

add_target_annotation <- function(top10_tbl) {
  if (is.null(top10_tbl) || nrow(top10_tbl) == 0) return(tibble())

  top10_tbl %>%
    left_join(
      PortalCompounds %>%
        select(
          CompoundID,
          GeneSymbolOfTargets,
          TargetOrMechanism
        ),
      by = c("ThirdDrugCompoundID" = "CompoundID")
    ) %>%
    mutate(TargetGene = map(GeneSymbolOfTargets, extract_targets)) %>%
    unnest(TargetGene, keep_empty = TRUE) %>%
    mutate(
      TargetGene = ifelse(is.na(TargetGene) | TargetGene == "", NA_character_, TargetGene)
    )
}

get_crispr_dependency <- function(model_id, target_gene) {

  if (!has_crispr || nrow(CRISPR_Dep) == 0) {
    return(tibble(
      CRISPR_Dependency = NA_real_,
      DependencyClass = "CRISPR file not loaded"
    ))
  }

  if (is.na(target_gene) || target_gene == "") {
    return(tibble(
      CRISPR_Dependency = NA_real_,
      DependencyClass = "No target annotation"
    ))
  }

  if (!target_gene %in% names(CRISPR_Dep)) {
    return(tibble(
      CRISPR_Dependency = NA_real_,
      DependencyClass = "Not available"
    ))
  }

  dep_value <- CRISPR_Dep %>%
    filter(ModelID == model_id) %>%
    pull(all_of(target_gene))

  dep_value <- suppressWarnings(as.numeric(dep_value))
  dep_value <- ifelse(length(dep_value) == 0, NA_real_, dep_value[1])

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
    return(tibble(
      Expression_log2TPM = NA_real_,
      ExpressionPercentile = NA_real_,
      ExpressionClass = "Expression file not loaded"
    ))
  }

  if (length(target_gene) == 0 || is.null(target_gene) ||
      is.na(target_gene[[1]]) || target_gene[[1]] == "") {
    return(tibble(
      Expression_log2TPM = NA_real_,
      ExpressionPercentile = NA_real_,
      ExpressionClass = "No target annotation"
    ))
  }

  target_gene <- target_gene[[1]]

  if (!target_gene %in% names(Expression_TPM)) {
    return(tibble(
      Expression_log2TPM = NA_real_,
      ExpressionPercentile = NA_real_,
      ExpressionClass = "Gene not available"
    ))
  }

  model_value <- Expression_TPM %>%
    filter(DepMapModelID == model_id) %>%
    pull(all_of(target_gene))

  model_value <- suppressWarnings(as.numeric(model_value))

  if (length(model_value) == 0) {
    return(tibble(
      Expression_log2TPM = NA_real_,
      ExpressionPercentile = NA_real_,
      ExpressionClass = "Model not available"
    ))
  }

  value <- model_value[[1]]

  if (is.na(value) || !is.finite(value)) {
    return(tibble(
      Expression_log2TPM = NA_real_,
      ExpressionPercentile = NA_real_,
      ExpressionClass = "Expression not available"
    ))
  }

  all_values <- suppressWarnings(as.numeric(Expression_TPM[[target_gene]]))
  all_values <- all_values[is.finite(all_values)]

  percentile <- if (length(all_values) > 0) {
    100 * mean(all_values <= value)
  } else {
    NA_real_
  }

  expr_class <- case_when(
    is.na(percentile) ~ "Not available",
    percentile >= 75 ~ "High expression",
    percentile >= 25 ~ "Moderate expression",
    TRUE ~ "Low expression"
  )

  tibble(
    Expression_log2TPM = value,
    ExpressionPercentile = percentile,
    ExpressionClass = expr_class
  )
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

  if (length(target_gene) == 0 || is.null(target_gene) ||
      is.na(target_gene[[1]]) || target_gene[[1]] == "") {
    return(tibble(
      MutationStatus = "No target annotation",
      MutationCount = NA_integer_,
      VariantClassification = NA_character_,
      ProteinChange = NA_character_
    ))
  }

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
    hits[[mutation_class_col]] %>%
      as.character() %>%
      na.omit() %>%
      unique() %>%
      paste(collapse = "; ")
  } else {
    NA_character_
  }

  protein_text <- if (!is.na(mutation_protein_col)) {
    hits[[mutation_protein_col]] %>%
      as.character() %>%
      na.omit() %>%
      unique() %>%
      paste(collapse = "; ")
  } else {
    NA_character_
  }

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
    return(list(
      functional = tibble(),
      molecular = tibble(),
      biomarker = tibble()
    ))
  }

  functional_tbl <- target_tbl %>%
    mutate(
      Functional_Profile = map2(ModelID, TargetGene, get_crispr_dependency)
    ) %>%
    unnest(Functional_Profile)

  molecular_tbl <- functional_tbl %>%
    mutate(
      Molecular_Profile = map2(ModelID, TargetGene, get_target_expression)
    ) %>%
    unnest(Molecular_Profile)

  biomarker_tbl <- molecular_tbl %>%
    mutate(
      Biomarker_Profile = map2(ModelID, TargetGene, get_target_mutation)
    ) %>%
    unnest(Biomarker_Profile)

  list(
    functional = functional_tbl,
    molecular = molecular_tbl,
    biomarker = biomarker_tbl
  )
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
          "Screen all GDSC2 compounds" = "all",
          "Select third drugs" = "selected"
        ),
        selected = "all"
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
          "Select one or multiple third drugs. Each selected drug is evaluated independently against the selected clinical reference backbone."
        )
      ),

      hr(),

      checkboxGroupInput(
        "doses",
        "Dose levels:",
        choices = c(0.001, 0.01, 0.1, 1),
        selected = c(0.001, 0.01, 0.1, 1)
      ),

      numericInput(
        "noise",
        "Noise SD:",
        value = 5,
        min = 0,
        max = 20,
        step = 1
      ),

      actionButton("run", "Run analysis", class = "btn-primary"),

      br(), br(),

      downloadButton("download_reference", "Download reference summary"),
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
        tabPanel("Third-drug ranking", DTOutput("ranking_table")),
        tabPanel("Top candidates per model", DTOutput("top10_table")),
        tabPanel("Bliss plot", plotOutput("top10_plot", height = "850px")),

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
          h4("Is dependency or Bliss improvement associated with a genomic alteration?"),
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
          mutate(SelectionStatus = "All available cervical GDSC2 compounds")
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

    validate(
      need(length(input$doses) > 0, "Please select at least one dose level.")
    )

    combo_row <- selected_combo()
    selected_drugs <- selected_third_drugs()

    if (input$third_drug_mode == "selected") {
      validate(
        need(nrow(selected_drugs) > 0, "Please select at least one third drug.")
      )
      selected_partner_ids <- unique(selected_drugs$CompoundID)
    } else {
      selected_partner_ids <- NULL
    }

    withProgress(
      message = "Running pharmacology and validation layers...",
      value = 0,
      {
        incProgress(0.25, detail = "Running Bliss analysis")

        out <- run_reference_pipeline(
          combo_name = combo_row$ClinicalCombination[[1]],
          therapeutic_class = combo_row$TherapeuticClass[[1]],
          drug_names = combo_row$DrugNames[[1]],
          compound_ids = combo_row$CompoundIDs[[1]],
          doses = as.numeric(input$doses),
          noise_sd = input$noise,
          selected_partner_ids = selected_partner_ids
        )

        incProgress(0.25, detail = "Adding CRISPR functional evidence")

        validation <- build_validation_tables(out$top10)

        out$third_drug_selection <- selected_drugs
        out$functional_tbl <- validation$functional

        incProgress(0.25, detail = "Adding target-expression evidence")

        out$molecular_tbl <- validation$molecular

        incProgress(0.15, detail = "Adding mutation biomarker context")

        out$biomarker_tbl <- validation$biomarker

        incProgress(0.10, detail = "Done")

        out
      }
    )
  })

  output$run_status <- renderText({
    if (is.null(results())) {
      return("Click 'Run analysis' to start.")
    }

    out <- results()

    screening_mode <- switch(
      input$third_drug_mode,
      all = "All cervical GDSC2 compounds",
      selected = "User-selected third drugs"
    )

    paste(
      paste0("Screening mode: ", screening_mode),
      paste0("Selected third drugs: ", nrow(out$third_drug_selection)),
      paste0("Eligible models: ", nrow(out$eligible_models)),
      paste0("Reference summary rows: ", nrow(out$reference_summary)),
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
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- tibble(Message = msg)
    }

    datatable(
      dat,
      filter = "top",
      options = list(pageLength = page_length, scrollX = TRUE)
    )
  }

  output$eligible_table <- renderDT({
    req(results())
    safe_dt(results()$eligible_models, "No eligible models.", 10)
  })

  output$reference_table <- renderDT({
    req(results())
    safe_dt(results()$reference_summary, "No reference summary available.", 10)
  })

  output$third_drug_selection_table <- renderDT({
    req(results())
    safe_dt(results()$third_drug_selection, "No third drugs selected.", 25)
  })

  output$ranking_table <- renderDT({
    req(results())
    safe_dt(results()$third_ranking, "No third-drug ranking available.", 25)
  })

  output$top10_table <- renderDT({
    req(results())
    safe_dt(results()$top10, "No candidate drugs available.", 25)
  })

  output$functional_table <- renderDT({
    req(results())
    safe_dt(
      results()$functional_tbl,
      "No functional CRISPR dependency evidence available.",
      25
    )
  })

  output$molecular_table <- renderDT({
    req(results())
    safe_dt(
      results()$molecular_tbl,
      "No molecular target-expression evidence available.",
      25
    )
  })

  output$biomarker_table <- renderDT({
    req(results())
    safe_dt(
      results()$biomarker_tbl,
      "No optional mutation biomarker evidence available.",
      25
    )
  })

  output$top10_plot <- renderPlot({
    req(results())

    plot_data <- results()$top10

    validate(
      need(nrow(plot_data) > 0, "No eligible third-drug results are available for plotting.")
    )

    reference_lines <- results()$reference_summary %>%
      distinct(ModelID, CellLineName, ReferenceMeanBliss)

    ggplot(
      plot_data,
      aes(
        x = reorder(ThirdDrug, MeanBliss),
        y = MeanBliss,
        fill = Interaction
      )
    ) +
      geom_col() +
      geom_hline(
        data = reference_lines,
        aes(yintercept = ReferenceMeanBliss),
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
          " + top third drugs"
        ),
        subtitle = "Dashed line = model-specific reference Mean Bliss",
        x = "Third drug",
        y = "Mean Bliss score",
        fill = "Interaction"
      )
  })

  output$functional_plot <- renderPlot({
    req(results())

    plot_data <- results()$functional_tbl %>%
      filter(
        !is.na(CRISPR_Dependency),
        DependencyClass %in% c(
          "Strong dependency", "Moderate dependency", "Weak dependency"
        )
      )

    validate(
      need(nrow(plot_data) > 0, "No CRISPR dependency data available for plotted targets.")
    )

    ggplot(
      plot_data,
      aes(
        x = reorder(TargetGene, CRISPR_Dependency),
        y = CRISPR_Dependency,
        fill = DependencyClass
      )
    ) +
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

    plot_data <- results()$molecular_tbl %>%
      filter(!is.na(Expression_log2TPM))

    validate(
      need(nrow(plot_data) > 0, "No expression data available for plotted targets.")
    )

    ggplot(
      plot_data,
      aes(
        x = reorder(TargetGene, Expression_log2TPM),
        y = Expression_log2TPM,
        fill = ExpressionClass
      )
    ) +
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

    plot_data <- results()$biomarker_tbl %>%
      count(CellLineName, MutationStatus, name = "n")

    validate(
      need(nrow(plot_data) > 0, "No mutation biomarker data available.")
    )

    ggplot(
      plot_data,
      aes(
        x = CellLineName,
        y = n,
        fill = MutationStatus
      )
    ) +
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

  output$download_reference <- downloadHandler(
    filename = function() paste0(input$combo, "_Reference_Summary.csv"),
    content = function(file) write.csv(results()$reference_summary, file, row.names = FALSE)
  )

  output$download_ranking <- downloadHandler(
    filename = function() paste0(input$combo, "_ThirdDrug_Ranking.csv"),
    content = function(file) write.csv(results()$third_ranking, file, row.names = FALSE)
  )

  output$download_top10 <- downloadHandler(
    filename = function() paste0(input$combo, "_Top_Candidates_Per_Model.csv"),
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

shinyApp(ui, server)
