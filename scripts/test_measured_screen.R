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
cervical_candidates <- intersect(Antiangiogenic_Compounds$CompoundID, Measured_Monotherapy$CompoundID)
antiangiogenic <- run_measured_screen(reference_ids, cervical_candidates, c(.03, .1))
unrestricted <- run_measured_screen(reference_ids, unique(Measured_Monotherapy$CompoundID),
                                    c(.03, .1))
stopifnot(nrow(antiangiogenic$ranking) > 0L, nrow(unrestricted$ranking) == 3062L)
cat("Measured-screen and Shiny server checks passed.\n")
example <- read_combination_benchmark(benchmark_example)
stopifnot(nrow(example) == 54L, dplyr::n_distinct(example$Model) == 6,
          all(table(example$Model) == 9L), all(example$BlissExpected >= example$HSAExpected),
          dplyr::n_distinct(Ovarian_Monotherapy$ModelID) == 41)
upload <- tempfile(fileext = ".csv")
writeLines(c(paste(benchmark_columns, collapse = ","), "Model1,run1,Drug1,Drug2,1,2,0.8,0.5,0.3"), upload)
shiny::testServer(combination_benchmark_server, {
  session$setInputs(benchmark_model = "OVCAR-4")
  stopifnot(nrow(selected_model()) == 9L)
  session$setInputs(benchmark_upload = list(datapath = upload), benchmark_model = "Model1")
  stopifnot(isTRUE(all.equal(selected_model()$BlissDeviationPP, 10)),
            isTRUE(all.equal(selected_model()$HSADeviationPP, 20)))
  stopifnot(!is.null(output$benchmark_plot))
})
shiny::testServer(server, {
  session$setInputs(cancer_type = "Ovary", combo = "Cisplatin_Paclitaxel",
                    candidate_mode = "all", doses = c("0.03","0.1"), run = 1)
  ovarian <- result()
  stopifnot(nrow(ovarian$ranking) > 0, all(ovarian$ranking$CancerType == "Ovary"))
  print(c(models = nrow(ovarian$reference), pairs = nrow(ovarian$ranking)))
})
