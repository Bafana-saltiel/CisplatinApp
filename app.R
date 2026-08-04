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
  "R/09_validation_agent.R",
  "R/10_cohort_validation.R",
  "R/07_ui.R",
  "R/08_server.R"
)

missing_modules <- module_files[!file.exists(module_files)]
if (length(missing_modules) > 0) {
  stop("Missing application modules: ", paste(missing_modules, collapse = ", "))
}

# Load every module into the same environment used to evaluate app.R.
# This is important on Posit Connect Cloud, where app.R may be evaluated in a
# managed application environment rather than directly in .GlobalEnv.
app_environment <- environment()
for (module_file in module_files) {
  sys.source(module_file, envir = app_environment)
}

required_app_objects <- c(
  "cross_scc_validation_ui",
  "cross_scc_validation_server",
  "ui",
  "server"
)
missing_app_objects <- required_app_objects[
  !vapply(
    required_app_objects,
    exists,
    logical(1),
    envir = app_environment,
    mode = "function",
    inherits = FALSE
  )
]
if (length(missing_app_objects) > 0) {
  stop(
    "Application modules did not define required functions: ",
    paste(missing_app_objects, collapse = ", ")
  )
}

shiny::shinyApp(ui = ui, server = server)
