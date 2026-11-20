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

rsconnect::writeManifest(appDir = ".")
message("manifest.json regenerated successfully.")
