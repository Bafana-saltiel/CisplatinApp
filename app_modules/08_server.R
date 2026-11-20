server <- function(input, output, session) {
  cross_scc_validation_server(input, output, session)
  measured_ids <- unique(Measured_Monotherapy$CompoundID)
  drug_choices <- Third_Drug_Choices[Third_Drug_Choices %in% measured_ids]
  updateSelectizeInput(session, "candidates", choices = drug_choices, server = TRUE)
  updateSelectInput(session, "mono_drug", choices = drug_choices, selected = "DPC-001793")

  annotate <- function(rows) {
    rows <- left_join(rows, select(Cervix_Model, ModelID, CellLineName), by = "ModelID")
    if ("CompoundID" %in% names(rows)) {
      rows <- left_join(rows, distinct(select(Available_Third_Drugs, CompoundID, CompoundName)),
                        by = "CompoundID")
    }
    rows
  }

  result <- eventReactive(input$run, {
    validate(need(length(input$doses), "Choose at least one dose."))
    regimen <- filter(Clinical_Cisplatin_Combinations_Final, ClinicalCombination == input$combo)
    req(nrow(regimen) == 1L)
    candidates <- switch(input$candidate_mode,
                         antiangiogenic = Cervical_Antiangiogenic_Compounds$CompoundID,
                         all = measured_ids,
                         selected = input$candidates)
    validate(need(length(candidates), "Select at least one candidate compound."))
    doses <- sort(as.numeric(input$doses))
    result <- run_measured_screen(regimen$CompoundIDs[[1]], candidates, doses)
    for (table in names(result)) {
      result[[table]] <- annotate(result[[table]]) %>%
        mutate(ReferenceRegimen = input$combo, CandidateGridUM = paste(doses, collapse = "|"))
    }
    result
  })

  output$ranking <- renderDT({
    rows <- result()$ranking
    validate(need(nrow(rows), "No model has measurements for the full reference and candidate dose grid. See Dose coverage."))
    datatable(rows, rownames = FALSE)
  })
  output$reference <- renderDT({
    rows <- result()$reference
    validate(need(nrow(rows), "No model has measurements for every component of this regimen."))
    datatable(rows, rownames = FALSE)
  })
  output$reference_doses <- renderDT(datatable(result()$reference_doses, rownames = FALSE))
  output$target_context <- renderDT({
    ranking <- result()$ranking %>%
      group_by(ModelID) %>%
      slice_head(n = 10L) %>%
      ungroup() %>%
      rename(ThirdDrugCompoundID = CompoundID)
    validate(need(nrow(ranking), "No ranked candidates have complete dose coverage."))
    context <- build_validation_tables(ranking)$biomarker
    datatable(context, rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$trials <- renderDT({
    req(input$combo)
    trials <- Trial_Arm_Matches %>% filter(Regimen == input$combo)
    validate(need(nrow(trials), "No same-arm registry match was identified for this preset."))
    datatable(trials, rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$coverage <- renderDT(datatable(result()$coverage, rownames = FALSE))

  monotherapy <- reactive({
    req(input$mono_drug)
    annotate(filter(Measured_Monotherapy, CompoundID == input$mono_drug))
  })
  output$measurements <- renderDT(datatable(monotherapy(), rownames = FALSE))
  output$measured_plot <- renderPlot({
    ggplot(monotherapy(), aes(DoseUM, 100 * ViabilityFraction, colour = CellLineName)) +
      geom_point() + scale_x_log10() + theme_minimal() +
      labs(x = "Measured concentration (µM)", y = "Vehicle-relative viability (%)", colour = "Model")
  })
  output$audit_summary <- renderDT({
    rows <- Monotherapy_Fit_Audit %>%
      pivot_longer(c(OriginalFit, SlopeReorientedFit), names_to = "Model", values_to = "Prediction") %>%
      group_by(Model) %>%
      summarise(N = n(), MAE = mean(abs(Prediction - ViabilityFraction)),
                RMSE = sqrt(mean((Prediction - ViabilityFraction)^2)),
                Bias = mean(Prediction - ViabilityFraction), .groups = "drop")
    datatable(rows, rownames = FALSE)
  })
  output$audit_table <- renderDT(datatable(Monotherapy_Fit_Audit, rownames = FALSE))
  output$download_audit <- downloadHandler(
    filename = function() "GDSC2MonotherapyFitAudit_Cervix.csv",
    content = function(path) write_csv(Monotherapy_Fit_Audit, path)
  )
  output$download_ranking <- downloadHandler(
    filename = function() "Cisplatin_BlissExpected_Increments.csv",
    content = function(path) write_csv(result()$ranking, path)
  )
}
