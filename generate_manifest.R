# Run this script locally from the application root before committing/publishing.
required_packages <- c(
  "shiny", "shinythemes", "DT", "dplyr", "stringr", "purrr",
  "tidyr", "ggplot2", "readr", "tibble", "igraph", "scales", "rsconnect"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# torch is optional. Install it locally only when you want Connect Cloud to try
# restoring and using the GCN engine. Otherwise the app uses graph propagation.
rsconnect::writeManifest(appDir = ".")
message("manifest.json regenerated successfully.")
