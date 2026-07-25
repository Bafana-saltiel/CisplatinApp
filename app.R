# Cisplatin Clinical Reference Combination Explorer
# Modular version with pharmacology, multi-synergy, validation, and graph ML.

module_files <- c(
  "R/00_packages.R",
  "R/01_data.R",
  "R/02_reference_regimens.R",
  "R/03_pharmacology.R",
  "R/04_compounds.R",
  "R/05_validation.R",
  "R/06_graph_ml.R",
  "R/07_ui.R",
  "R/08_server.R"
)

missing_modules <- module_files[!file.exists(module_files)]
if (length(missing_modules) > 0) {
  stop("Missing application modules: ", paste(missing_modules, collapse = ", "))
}

invisible(lapply(module_files, source, local = FALSE))
shiny::shinyApp(ui = ui, server = server)
