# Run from the application root.
app <- source("app.R")$value
fixture <- tibble(
  ModelID = "test",
  CompoundID = rep(c("a", "b", "c", "outside"), each = 2),
  DoseUM = c(1, 2, 1, 2, 1, 2, 2, 3),
  ViabilityFraction = c(.5, .5, .6, .6, .4, .4, .1, .1)
)
fixture_result <- run_measured_screen(c("a", "b"), c("c", "outside"), 1, fixture)
stopifnot(
  identical(fixture_result$ranking$CompoundID, "c"),
  isTRUE(all.equal(fixture_result$reference$ReferenceBlissExpectedInhibitionPct, 70)),
  isTRUE(all.equal(fixture_result$ranking$ExpandedBlissExpectedInhibitionPct, 88)),
  isTRUE(all.equal(fixture_result$ranking$BlissExpectedIncrementPP, 18)),
  !fixture_result$coverage$GridComplete[fixture_result$coverage$CompoundID == "outside"],
  nrow(run_measured_screen(c("a", "missing"), "c", 1, fixture)$ranking) == 0L
)
stopifnot(isTRUE(all.equal(
  estimate_at_doses(tibble(DoseUM = c(1, 100), ViabilityFraction = c(.8, .4)), c(.1, 10, 100)),
  c(NA_real_, .6, .4)
)))

shiny::testServer(server, {
  session$setInputs(combo = "Cisplatin_Paclitaxel", candidate_mode = "selected",
                    candidates = "DPC-001165", doses = c("0.03", "0.1"),
                    run = 1)
  data <- result()
  stopifnot(
    nrow(data$reference) == 13L,
    nrow(data$ranking) > 0L,
    all(data$ranking$CompoundID == "DPC-001165"),
    any(data$ranking$BlissExpectedIncrementPP > 0),
    nrow(data$coverage) >= nrow(data$ranking),
    nrow(monotherapy()) > 0L
  )
  invisible(output$ranking)
  invisible(output$reference)
  stopifnot(all(monotherapy()$CompoundID == "DPC-001165"),
            nchar(output$measured_plot$src) > 1000)
  invisible(output$measurements)
  invisible(output$audit_summary)
  invisible(output$trials)
  invisible(output$target_context)
  download <- output$download_ranking
  stopifnot(file.exists(download), nrow(read_csv(download, show_col_types = FALSE)) == nrow(data$ranking))
})
reference_ids <- c("DPC-001793", "DPC-004880")
antiangiogenic <- run_measured_screen(reference_ids, Cervical_Antiangiogenic_Compounds$CompoundID,
                                      c(.03, .1))
unrestricted <- run_measured_screen(reference_ids, unique(Measured_Monotherapy$CompoundID),
                                    c(.03, .1))
stopifnot(nrow(antiangiogenic$ranking) > 0L, nrow(unrestricted$ranking) == 3062L)
cat("Measured-screen and Shiny server checks passed.\n")
