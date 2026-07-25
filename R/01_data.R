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

