benchmark_example <- file.path(DATA_DIR, "ALMANAC_ovarian_cisplatin_paclitaxel.csv")
benchmark_columns <- c("Model", "Experiment", "DrugA", "DrugB", "DoseA_uM", "DoseB_uM",
                       "ViabilityA", "ViabilityB", "CombinationViability")

read_combination_benchmark <- function(path) {
  rows <- readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(
    .default = "c", DoseA_uM = "d", DoseB_uM = "d", ViabilityA = "d",
    ViabilityB = "d", CombinationViability = "d"
  ))
  if (!all(benchmark_columns %in% names(rows)) || !nrow(rows) ||
      nrow(readr::problems(rows)) || anyNA(rows[benchmark_columns])) {
    stop("Enter one complete row per measured dose pair using the example CSV columns.")
  }
  numeric_columns <- benchmark_columns[5:9]
  if (any(!is.finite(as.matrix(rows[numeric_columns]))) ||
      any(rows$DoseA_uM <= 0 | rows$DoseB_uM <= 0) ||
      any(rows$ViabilityA < 0 | rows$ViabilityB < 0 | rows$CombinationViability < 0) ||
      any(tolower(rows$DrugA) == tolower(rows$DrugB))) {
    stop("Use positive doses in µM, non-negative viability fractions and distinct drugs.")
  }
  rows %>%
    mutate(
      ObservedInhibition = 100 * (1 - CombinationViability),
      BlissExpected = 100 * (1 - pmin(ViabilityA, 1) * pmin(ViabilityB, 1)),
      HSAExpected = 100 * (1 - pmin(ViabilityA, ViabilityB, 1)),
      BlissDeviationPP = ObservedInhibition - BlissExpected,
      HSADeviationPP = ObservedInhibition - HSAExpected
    )
}

combination_benchmark_ui <- function() {
  tabPanel(
    "Experimental combinations",
    p("The built-in example uses NCI ALMANAC ovarian cisplatin–paclitaxel measurements. Upload a CSV to analyse your own two-drug experiment."),
    p("Enter one row per tested dose pair, with both matching single-agent measurements and the observed combination response. Use vehicle-relative fractions (1 = vehicle control) and doses in µM. Single-agent fractions above 1 are bounded at 1 for both expectations. Observed combination values are retained as entered."),
    downloadButton("benchmark_example", "Download example CSV"),
    fileInput("benchmark_upload", "Upload experimental dose pairs", accept = ".csv"),
    uiOutput("benchmark_model_selector"),
    plotOutput("benchmark_plot"),
    DTOutput("benchmark_results"),
    downloadButton("benchmark_download", "Download computed comparisons")
  )
}

combination_benchmark_server <- function(input, output, session) {
  comparisons <- reactive({
    path <- if (is.null(input$benchmark_upload)) benchmark_example else input$benchmark_upload$datapath
    rows <- tryCatch(read_combination_benchmark(path), error = function(error) error)
    validate(need(!inherits(rows, "error"), if (inherits(rows, "error")) conditionMessage(rows)))
    rows
  })
  output$benchmark_model_selector <- renderUI({
    selectInput("benchmark_model", "Experimental model", choices = sort(unique(comparisons()$Model)))
  })
  selected_model <- reactive({
    req(input$benchmark_model)
    rows <- comparisons() %>% filter(Model == input$benchmark_model)
    req(nrow(rows))
    rows
  })
  output$benchmark_plot <- renderPlot({
    rows <- selected_model() %>%
      pivot_longer(c(BlissExpected, HSAExpected), names_to = "Benchmark", values_to = "Expected")
    ggplot(rows, aes(Expected, ObservedInhibition, colour = Benchmark)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(size = 2.5, alpha = 0.75) +
      coord_equal() +
      scale_colour_discrete(labels = c(BlissExpected = "Bliss independence", HSAExpected = "HSA")) +
      labs(x = "Expected inhibition (%)", y = "Observed inhibition (%)", colour = NULL) +
      theme_bw()
  })
  output$benchmark_results <- renderDT({
    datatable(selected_model(), rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$benchmark_example <- downloadHandler(
    filename = function() basename(benchmark_example),
    content = function(path) file.copy(benchmark_example, path)
  )
  output$benchmark_download <- downloadHandler(
    filename = function() "experimental_combination_comparisons.csv",
    content = function(path) write_csv(comparisons(), path)
  )
}
