# ============================================================
# Shiny App: Cisplatin Clinical Reference Combination Explorer
# Updated with monotherapy response, selectable synergy models, anti-angiogenic screening, graph neural network prioritisation, and validation layers
# ============================================================

required_packages <- c(
  "shiny", "shinythemes", "DT", "dplyr", "stringr",
  "purrr", "tidyr", "ggplot2", "readr", "tibble", "igraph", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Install required packages before starting: ", paste(missing_packages, collapse = ", "))
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
  library(igraph)
  library(scales)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
