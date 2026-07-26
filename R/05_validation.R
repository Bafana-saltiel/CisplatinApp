# ============================================================
# 6. TARGET-LEVEL VALIDATION MODULES
# ============================================================

extract_target_genes <- function(target_string) {
  if (is.null(target_string) || length(target_string) == 0) {
    return(character(0))
  }

  value <- as.character(target_string[[1]])
  if (is.na(value) || !nzchar(trimws(value))) {
    return(character(0))
  }

  genes <- unlist(strsplit(gsub(",", ";", value, fixed = TRUE), ";", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  unique(genes)
}

add_target_annotation <- function(top10_tbl) {
  if (is.null(top10_tbl) || nrow(top10_tbl) == 0) return(tibble())

  top10_tbl %>%
    left_join(
      PortalCompounds %>% select(CompoundID, GeneSymbolOfTargets, TargetOrMechanism),
      by = c("ThirdDrugCompoundID" = "CompoundID")
    ) %>%
    mutate(TargetGene = map(GeneSymbolOfTargets, extract_target_genes)) %>%
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
