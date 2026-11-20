# Run from the application root with Rscript scripts/reproduce_screen.R OUTPUT_DIR.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide an output directory.")
dir.create(args[[1]], recursive = TRUE, showWarnings = FALSE)
source("app.R")

doses <- c(0.03, 0.1)
screen <- run_measured_screen(
  c("DPC-001793", "DPC-004880"), unique(Measured_Monotherapy$CompoundID), doses
)
candidate_ids <- c("DPC-000247", "DPC-001165", "DPC-000258", "DPC-002313",
                   "DPC-002900", "DPC-006002", "DPC-004707")
example <- screen$ranking %>%
  filter(CompoundID %in% candidate_ids) %>%
  left_join(select(Cervix_Model, ModelID, CellLineName), by = "ModelID") %>%
  left_join(select(Available_Third_Drugs, CompoundID, CompoundName), by = "CompoundID")
screen$seven_candidates <- example
for (name in names(screen)) {
  write_csv(screen[[name]], file.path(args[[1]], paste0(name, ".csv")))
}
writeLines(capture.output(sessionInfo()), file.path(args[[1]], "R-session-info.txt"))
