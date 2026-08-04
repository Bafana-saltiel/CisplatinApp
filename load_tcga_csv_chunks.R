# Load complete TCGA expression matrices stored as GitHub-safe CSV.gz chunks.

load_tcga_csv_chunks <- function(project, data_dir = "data/TCGA data", targets = NULL) {
  project <- toupper(project)
  project_dir <- file.path(data_dir, paste0("TCGA_", project))

  chunk_files <- sort(list.files(
    project_dir,
    pattern = paste0("^TCGA_", project, "_expression_part_[0-9]+[ab]?\\.csv\\.gz$"),
    full.names = TRUE
  ))

  if (length(chunk_files) == 0) {
    stop("No expression chunks found for TCGA-", project, " in ", project_dir)
  }

  expression_parts <- lapply(chunk_files, function(path) {
    part <- readr::read_csv(path, show_col_types = FALSE)
    if (!is.null(targets)) {
      part <- dplyr::filter(part, GeneSymbol %in% targets)
    }
    part
  })

  expression <- dplyr::bind_rows(expression_parts)
  if (is.null(targets) && nrow(expression) != 60660) {
    stop("Expected 60,660 genes for TCGA-", project, "; found ", nrow(expression))
  }

  metadata_file <- file.path(project_dir, paste0("TCGA_", project, "_sample_metadata.csv"))
  if (!file.exists(metadata_file)) {
    stop("Missing metadata file: ", metadata_file)
  }
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)

  expression_patients <- names(expression)[-(1:2)]
  if (!identical(expression_patients, metadata$PatientID)) {
    stop("Expression columns and metadata PatientID order do not match for TCGA-", project)
  }

  list(
    project = paste0("TCGA-", project),
    expression = expression,
    metadata = metadata,
    expression_unit = "log2(FPKM-UQ + 1)",
    gene_annotation = "GENCODE v36"
  )
}

# Examples:
# cesc_complete <- load_tcga_csv_chunks("CESC")
# cesc_targets <- load_tcga_csv_chunks("CESC", targets = CROSS_SCC_TARGETS)
