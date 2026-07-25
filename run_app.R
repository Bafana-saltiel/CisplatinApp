# Local launcher. Install dependencies before running this file.
required_packages <- c(
  "shiny", "shinythemes", "DT", "dplyr", "stringr", "purrr",
  "tidyr", "ggplot2", "readr", "tibble", "igraph", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the missing package(s) before starting the app: ",
    paste(missing_packages, collapse = ", ")
  )
}
shiny::runApp(".")
