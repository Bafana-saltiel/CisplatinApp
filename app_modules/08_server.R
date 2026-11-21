server <- function(input, output, session) {
  measured_ids <- unique(Measured_Monotherapy$CompoundID)
  drug_choices <- Third_Drug_Choices[Third_Drug_Choices %in% measured_ids]
  updateSelectizeInput(session, "candidates", choices = drug_choices, server = TRUE)

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
    result <- run_measured_screen(regimen$CompoundIDs[[1]], candidates, doses,
                                  include_hsa = isTRUE(input$include_hsa),
                                  include_loewe = isTRUE(input$include_loewe),
                                  include_zip = isTRUE(input$include_zip))
    for (table in names(result)) {
      result[[table]] <- annotate(result[[table]]) %>%
        mutate(ReferenceRegimen = input$combo, CandidateGridUM = paste(doses, collapse = "|"))
    }
    result
  })

  output$hsa_note <- renderUI({
    metrics <- names(result()$ranking)
    tagList(
      if ("HSAExpectedIncrementPP" %in% metrics)
        p("HSA averages the strongest single-agent inhibition at each dose combination. Candidate ranking uses Bliss."),
      if ("LoeweExpectedIncrementPP" %in% metrics)
        p("Loewe increments use matched dose tuples where both expectations are defined within the fitted monotherapy effect ranges. Coverage is the fraction of tuples used."),
      if ("ZIPExpectedIncrementPP" %in% metrics)
        p("ZIP expectation combines fitted single-agent responses under zero potency shift. ZIP synergy requires measured combination responses.")
    )
  })
  output$ranking <- renderDT({
    rows <- result()$ranking
    validate(need(nrow(rows), "No model has measurements for the full reference and candidate dose grid. See Dose coverage."))
    datatable(rows, rownames = FALSE)
  })
  top_candidates <- reactive({
    rows <- result()$ranking
    validate(need(nrow(rows), "No eligible candidates on the selected dose grid."))
    rows %>%
      arrange(CellLineName, ModelID, desc(BlissExpectedIncrementPP), CompoundID) %>%
      group_by(ModelID) %>%
      slice_head(n = 10L) %>%
      ungroup()
  })

  ranked_target_context <- reactive({
    ranking <- top_candidates() %>%
      group_by(ModelID) %>%
      mutate(CandidateRank = row_number()) %>%
      ungroup() %>%
      rename(ThirdDrugCompoundID = CompoundID)
    build_validation_tables(ranking)$biomarker
  })
  candidate_cohort_server(input, output, session, ranked_target_context)
  cross_scc_validation_server(input, output, session, ranked_target_context)

  output$top_candidates_plot <- renderPlot({
    rows <- top_candidates() %>%
      mutate(
        Panel = if_else(is.na(CellLineName), ModelID, CellLineName),
        Tick = factor(rev(seq_len(n()))),
        Comparison = if_else(BlissExpectedIncrementPP > 1e-8,
                             "Higher expected inhibition", "No expected increase")
      )
    labels <- setNames(
      if_else(is.na(rows$CompoundName), rows$CompoundID, rows$CompoundName),
      as.character(rows$Tick)
    )
    reference <- distinct(rows, Panel, ReferenceBlissExpectedInhibitionPct)

    ggplot(rows, aes(ExpandedBlissExpectedInhibitionPct, Tick, fill = Comparison)) +
      geom_col(width = 0.85) +
      geom_vline(data = reference, aes(xintercept = ReferenceBlissExpectedInhibitionPct),
                 linetype = "dashed", linewidth = 0.55) +
      facet_wrap(vars(Panel), ncol = 4, scales = "free_y") +
      scale_y_discrete(labels = labels) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
      scale_fill_manual(values = c("Higher expected inhibition" = "#009D46",
                                   "No expected increase" = "#638DD5")) +
      labs(
        title = paste(gsub("_", " + ", rows$ReferenceRegimen[[1]]), "with candidate additions"),
        subtitle = paste0("Dashed line shows the reference expectation for each model. Candidate grid (µM) ",
                          gsub("|", ", ", rows$CandidateGridUM[[1]], fixed = TRUE)),
        caption = "Bliss independence expectations calculated from measured monotherapy.",
        x = "Bliss-expected inhibition (%)", y = "Candidate compound", fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
  }, height = function() {
    max(360, 110 + 225 * ceiling(n_distinct(top_candidates()$ModelID) / 4))
  })

  output$reference <- renderDT({
    rows <- result()$reference
    validate(need(nrow(rows), "No model has measurements for every component of this regimen."))
    datatable(rows, rownames = FALSE)
  })
  output$reference_doses <- renderDT(datatable(result()$reference_doses, rownames = FALSE))
  target_gene_rows <- function(score, evidence_class) {
    rows <- ranked_target_context()
    validate(need(nrow(rows), "No target context is available for this screen."))
    rows <- rows %>%
      filter(!is.na(TargetGene), is.finite(.data[[score]])) %>%
      distinct(ModelID, CellLineName, TargetGene, .data[[score]], .data[[evidence_class]]) %>%
      mutate(Panel = coalesce(CellLineName, ModelID)) %>%
      arrange(Panel, desc(.data[[score]]), TargetGene) %>%
      mutate(Tick = factor(rev(row_number())))
    validate(need(nrow(rows), "No measurements are available for the ranked candidate targets."))
    rows
  }

  render_target_gene_plot <- function(score, evidence_class, title, subtitle, axis, palette,
                                      ncol, thresholds = numeric()) {
    renderPlot({
      rows <- target_gene_rows(score, evidence_class)
      plot <- ggplot(rows, aes(.data[[score]], Tick, fill = .data[[evidence_class]])) +
        geom_col(orientation = "y") +
        facet_wrap(vars(Panel), ncol = ncol, scales = "free_y") +
        scale_y_discrete(labels = setNames(rows$TargetGene, as.character(rows$Tick))) +
        scale_fill_manual(values = palette) +
        labs(title = title, subtitle = subtitle, x = axis, y = "Target gene",
             fill = if (length(thresholds)) "Dependency class" else "Expression class") +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 8), legend.position = "bottom")
      if (length(thresholds)) {
        plot <- plot +
          geom_vline(xintercept = thresholds[1], linetype = "dashed") +
          geom_vline(xintercept = thresholds[2], linetype = "dotted") +
          scale_x_continuous(breaks = seq(0, 1, .25)) +
          coord_cartesian(xlim = c(0, 1))
      }
      plot
    }, height = function() {
      counts <- as.integer(table(target_gene_rows(score, evidence_class)$Panel))
      facet_rows <- split(counts, ceiling(seq_along(counts) / ncol))
      110 + sum(vapply(facet_rows, function(x) max(150, 45 + 14 * max(x)), numeric(1)))
    })
  }

  output$target_crispr_plot <- render_target_gene_plot(
    "CRISPR_Dependency", "DependencyClass",
    "Functional evidence: CRISPR dependency of third-drug targets",
    "Dashed: moderate cutoff (0.50). Dotted: strong cutoff (0.80).",
    "CRISPR dependency probability",
    c("Moderate dependency" = "#F8766D", "Strong dependency" = "#00BA38",
      "Weak dependency" = "#619CFF", "Not dependent" = "#999999"),
    ncol = 3, thresholds = c(.5, .8)
  )

  output$target_expression_plot <- render_target_gene_plot(
    "Expression_log2TPM", "ExpressionClass",
    "Molecular evidence: target expression in matched cervical models",
    "Expression classes use target-specific percentiles.",
    "log2(TPM + 1)",
    c("High expression" = "#F8766D", "Low expression" = "#00BA38",
      "Moderate expression" = "#619CFF", "Not available" = "#999999"),
    ncol = 4
  )

  output$target_mutation_plot <- renderPlot({
    rows <- ranked_target_context()
    validate(need(nrow(rows), "No target context is available for this screen."))
    counts <- rows %>%
      distinct(ModelID, CellLineName, ThirdDrugCompoundID, TargetGene, MutationStatus) %>%
      mutate(
        Panel = coalesce(CellLineName, ModelID),
        Status = case_when(
          MutationStatus == "Mutated" ~ "Mutated",
          MutationStatus == "No target annotation" ~ "No target annotation",
          MutationStatus == "No recorded mutation" ~ "No recorded mutation",
          TRUE ~ "Mutation data unavailable"
        ),
        Status = factor(Status, levels = c("Mutated", "No target annotation",
                                           "Mutation data unavailable", "No recorded mutation"))
      ) %>%
      count(Panel, Status)
    ggplot(counts, aes(n, Panel, fill = Status)) +
      geom_col(orientation = "y") +
      scale_fill_manual(values = c("Mutated" = "#F8766D", "No target annotation" = "#00BA38",
                                   "No recorded mutation" = "#619CFF",
                                   "Mutation data unavailable" = "#999999")) +
      labs(title = "Somatic mutation context for candidate drug targets",
           x = "Drug–target associations or unannotated candidates",
           y = "Cervical cancer model", fill = "Mutation status") +
      theme_bw(base_size = 11)
  }, height = function() {
    max(350, 110 + 34 * n_distinct(ranked_target_context()$ModelID))
  })

  output$target_context <- renderDT({
    datatable(ranked_target_context(), rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$trials <- renderDT({
    req(input$combo)
    trials <- Trial_Arm_Matches %>% filter(Regimen == input$combo)
    validate(need(nrow(trials), "No same-arm registry match was identified for this preset."))
    datatable(trials, rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$coverage <- renderDT(datatable(result()$coverage, rownames = FALSE))

  top_monotherapy <- reactive({
    validate(need(input$run > 0, "Calculate expectations to see the top monotherapy responses."))
    rows <- result()$ranking
    validate(need(nrow(rows), "No eligible monotherapy responses on the selected dose grid."))
    rows %>%
      arrange(CellLineName, ModelID, desc(MeanSingleAgentInhibitionPct), CompoundID) %>%
      group_by(ModelID) %>%
      slice_head(n = 10L) %>%
      ungroup() %>%
      mutate(Panel = coalesce(CellLineName, ModelID),
             Tick = factor(rev(seq_len(n()))),
             Label = coalesce(CompoundName, CompoundID),
             MonotherapyClass = case_when(
               MeanSingleAgentInhibitionPct >= 75 ~ "High intrinsic response",
               MeanSingleAgentInhibitionPct >= 40 ~ "Moderate intrinsic response",
               MeanSingleAgentInhibitionPct > 0 ~ "Low intrinsic response",
               TRUE ~ "No positive inhibition"
             ))
  })

  output$top_monotherapy_plot <- renderPlot({
    rows <- top_monotherapy()
    labels <- setNames(rows$Label, as.character(rows$Tick))
    ggplot(rows, aes(MeanSingleAgentInhibitionPct, Tick, fill = MonotherapyClass)) +
      geom_col(data = filter(rows, MeanSingleAgentInhibitionPct > 0),
               width = 0.9, orientation = "y", show.legend = TRUE) +
      geom_point(data = filter(rows, MonotherapyClass == "No positive inhibition"),
                 aes(x = 0), shape = 21, colour = "#C77CFF", size = 2) +
      facet_wrap(vars(Panel), ncol = 4, scales = "free_y") +
      scale_y_discrete(labels = labels) +
      scale_x_continuous(breaks = seq(0, 100, 25)) +
      coord_cartesian(xlim = c(0, 100)) +
      scale_fill_manual(values = c("High intrinsic response" = "#F8766D",
                                   "Low intrinsic response" = "#7CAE00",
                                   "Moderate intrinsic response" = "#00BFC4",
                                   "No positive inhibition" = "#C77CFF")) +
      labs(title = "Top monotherapy responses in matched cervical models",
           subtitle = paste0("Measured GDSC2 viability summarised on the candidate grid (µM) ",
                             gsub("|", ", ", rows$CandidateGridUM[[1]], fixed = TRUE)),
           x = "Mean single-agent inhibition (%)", y = "Candidate compound",
           fill = "Monotherapy class",
           caption = "High: ≥75%. Moderate: 40 to <75%. Low: >0 to <40%. Purple dot: ≤0%.") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 8), legend.position = "bottom")
  }, height = function() {
    max(360, 110 + 245 * ceiling(n_distinct(top_monotherapy()$ModelID) / 4))
  })

  monotherapy <- reactive({
    pairs <- top_monotherapy() %>%
      select(ModelID, CompoundID, Panel, Label)
    Measured_Monotherapy %>%
      inner_join(pairs, by = c("ModelID", "CompoundID")) %>%
      mutate(InhibitionPct = 100 * (1 - ViabilityFraction)) %>%
      arrange(Panel, Label, CompoundID, DoseUM) %>%
      select(CellLineName = Panel, CompoundName = Label, everything())
  })
  output$measurements <- renderDT(datatable(monotherapy(), rownames = FALSE,
                                            options = list(scrollX = TRUE)))
  output$measured_plot <- renderPlot({
    rows <- monotherapy()
    line_rows <- rows %>% group_by(ModelID, CompoundID) %>% filter(n() > 1L) %>% ungroup()
    ggplot(rows, aes(DoseUM, InhibitionPct, colour = CompoundName,
                     group = interaction(ModelID, CompoundID))) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_line(data = line_rows, linewidth = 0.55) +
      geom_point(size = 1.4) +
      facet_wrap(vars(CellLineName), ncol = 4) +
      scale_x_log10() +
      labs(title = "Third-drug monotherapy dose–response profiles",
           subtitle = "Top ten monotherapy candidates per model at all available measured doses",
           x = "Measured concentration (µM)", y = "Measured monotherapy inhibition (%)",
           colour = "Third drug") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  }, height = function() {
    max(360, 170 + 230 * ceiling(n_distinct(monotherapy()$ModelID) / 4))
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
