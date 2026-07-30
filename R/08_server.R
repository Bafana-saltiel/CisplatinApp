# ============================================================
# 8. SHINY SERVER
# ============================================================

server <- function(input, output, session) {
  selected_combo <- reactive({
    Clinical_Cisplatin_Combinations_Final %>%
      filter(ClinicalCombination == input$combo)
  })

  selected_third_drugs <- reactive({
    if (input$third_drug_mode == "all") {
      return(
        Available_Third_Drugs %>%
          mutate(
            SelectionStatus =
              "All available cervical GDSC2 compounds"
          )
      )
    }

    if (input$third_drug_mode == "antiangiogenic") {
      return(
        Cervical_Antiangiogenic_Compounds %>%
          mutate(
            SelectionStatus =
              "Automatically selected anti-angiogenic cervical compound"
          )
      )
    }

    if (input$third_drug_mode == "selected") {
      req(input$selected_third_drugs)

      return(
        Available_Third_Drugs %>%
          filter(CompoundID %in% input$selected_third_drugs) %>%
          mutate(SelectionStatus = "User selected")
      )
    }

    tibble()
  })

  results <- eventReactive(input$run, {
    validate(need(length(input$doses) > 0, "Please select at least one dose level."))

    combo_row <- selected_combo()
    selected_drugs <- selected_third_drugs()

    if (input$third_drug_mode %in% c("antiangiogenic", "selected")) {
      empty_message <- if (input$third_drug_mode == "antiangiogenic") {
        paste(
          "No anti-angiogenic compounds with cervical GDSC2 response",
          "data were identified. Check the compound annotations and files."
        )
      } else {
        "Please select at least one third drug."
      }

      validate(need(nrow(selected_drugs) > 0, empty_message))
      selected_partner_ids <- unique(selected_drugs$CompoundID)
    } else {
      selected_partner_ids <- NULL
    }

    withProgress(message = "Running pharmacology and validation layers...", value = 0, {
      incProgress(0.20, detail = paste("Running reference-regimen", synergy_method_label(input$synergy_method), "analysis"))

      out <- run_reference_pipeline(
        combo_name = combo_row$ClinicalCombination[[1]],
        therapeutic_class = combo_row$TherapeuticClass[[1]],
        drug_names = combo_row$DrugNames[[1]],
        compound_ids = combo_row$CompoundIDs[[1]],
        doses = as.numeric(input$doses),
        noise_sd = input$noise,
        synergy_method = input$synergy_method,
        selected_partner_ids = selected_partner_ids
      )

      incProgress(0.20, detail = "Evaluating monotherapy response")
      out$third_drug_selection <- selected_drugs

      incProgress(0.20, detail = "Building target-level validation tables")
      validation <- build_validation_tables(out$top10)

      incProgress(0.15, detail = "Adding CRISPR functional evidence")
      out$functional_tbl <- validation$functional

      incProgress(0.15, detail = "Adding target-expression evidence")
      out$molecular_tbl <- validation$molecular

      incProgress(0.10, detail = "Adding mutation biomarker context")
      out$biomarker_tbl <- validation$biomarker

      if (isTRUE(input$enable_gnn)) {
        incProgress(0.10, detail = "Constructing drug graphs and training graph model")
        graph_out <- run_graph_learning(
          third_ranking = out$third_ranking,
          synergy_method = input$synergy_method,
          edge_threshold = input$graph_edge_threshold,
          epochs = input$gnn_epochs,
          hidden_dim = input$gnn_hidden_dim,
          engine = input$graph_engine
        )
        out$graph_ranking <- graph_out$ranking
        out$graph_diagnostics <- graph_out$diagnostics
        out$graph_edges <- graph_out$edges
        out$graph_features <- graph_out$features

        if (isTRUE(input$enable_validation_agent)) {
          incProgress(0.10, detail = "Running Validation Agent and graph-model benchmarking")
          validation_agent <- run_validation_agent(
            graph_features = graph_out$features,
            synergy_method = input$synergy_method,
            edge_threshold = input$graph_edge_threshold,
            requested_folds = input$validation_folds,
            decision_threshold = input$validation_threshold,
            epochs = input$gnn_epochs,
            hidden_dim = input$gnn_hidden_dim,
            engine = input$graph_engine
          )
          out$validation_predictions <- validation_agent$predictions
          out$validation_fold_metrics <- validation_agent$fold_metrics
          out$validation_summary <- validation_agent$summary
          out$validation_confusion <- validation_agent$confusion
          out$validation_roc <- validation_agent$roc
          out$validation_decision_curve <- validation_agent$decision_curve
          out$validation_clustering <- validation_agent$clustering
        } else {
          out$validation_predictions <- tibble()
          out$validation_fold_metrics <- tibble()
          out$validation_summary <- tibble()
          out$validation_confusion <- tibble()
          out$validation_roc <- tibble()
          out$validation_decision_curve <- tibble()
          out$validation_clustering <- tibble()
        }
      } else {
        out$graph_ranking <- tibble()
        out$graph_diagnostics <- tibble()
        out$graph_edges <- tibble()
        out$graph_features <- tibble()
        out$validation_predictions <- tibble()
        out$validation_fold_metrics <- tibble()
        out$validation_summary <- tibble()
        out$validation_confusion <- tibble()
        out$validation_roc <- tibble()
        out$validation_decision_curve <- tibble()
        out$validation_clustering <- tibble()
      }

      out
    })
  })

  observeEvent(results(), {
    out <- results()

    if (
      is.null(out$graph_ranking) ||
      nrow(out$graph_ranking) == 0 ||
      !all(c("ModelID", "CellLineName") %in% names(out$graph_ranking))
    ) {
      updateSelectInput(
        session,
        "graph_model_view",
        choices = character(0),
        selected = character(0)
      )
      return(invisible(NULL))
    }

    graph_models <- out$graph_ranking %>%
      distinct(ModelID, CellLineName) %>%
      arrange(CellLineName)

    choices <- setNames(graph_models$ModelID, graph_models$CellLineName)

    updateSelectInput(
      session,
      "graph_model_view",
      choices = choices,
      selected = unname(choices[[1]])
    )
  }, ignoreInit = TRUE)

  output$run_status <- renderText({
    if (is.null(results())) return("Click 'Run analysis' to start.")

    out <- results()
    screening_mode <- switch(
      input$third_drug_mode,
      all = "All cervical GDSC2 compounds",
      antiangiogenic = paste0(
        "All anti-angiogenic cervical compounds (n = ",
        nrow(Cervical_Antiangiogenic_Compounds),
        ")"
      ),
      selected = "User-selected third drugs"
    )

    paste(
      paste0("Screening mode: ", screening_mode),
      paste0("Synergy model: ", synergy_method_label(input$synergy_method)),
      paste0("Selected third drugs: ", nrow(out$third_drug_selection)),
      paste0("Eligible models: ", nrow(out$eligible_models)),
      paste0("Reference summary rows: ", nrow(out$reference_summary)),
      paste0("Monotherapy summary rows: ", nrow(out$monotherapy_summary)),
      paste0("Monotherapy dose-response rows: ", nrow(out$monotherapy_dose_response)),
      paste0("Third-drug ranking rows: ", nrow(out$third_ranking)),
      paste0("Top candidate rows: ", nrow(out$top10)),
      paste0("Functional CRISPR rows: ", nrow(out$functional_tbl)),
      paste0("Molecular expression rows: ", nrow(out$molecular_tbl)),
      paste0("Optional biomarker rows: ", nrow(out$biomarker_tbl)),
      paste0("Graph learning enabled: ", isTRUE(input$enable_gnn)),
      paste0("Graph ranking rows: ", nrow(out$graph_ranking)),
      paste0("Graph engine requested: ", ifelse(isTRUE(input$enable_gnn), input$graph_engine, "not run")),
      paste0("Validation Agent enabled: ", isTRUE(input$enable_validation_agent)),
      paste0("Cross-validation prediction rows: ", nrow(out$validation_predictions)),
      paste0("Validation summary rows: ", nrow(out$validation_summary)),
      paste0("R torch available: ", has_torch),
      paste0("CRISPR loaded: ", has_crispr),
      paste0("Expression loaded: ", has_expression),
      paste0("Mutation MAF loaded: ", has_mutation),
      sep = "\n"
    )
  })

  safe_dt <- function(dat, msg, page_length = 25) {
    if (is.null(dat) || nrow(dat) == 0) dat <- tibble(Message = msg)
    datatable(dat, filter = "top", options = list(pageLength = page_length, scrollX = TRUE))
  }

  output$eligible_table <- renderDT({ req(results()); safe_dt(results()$eligible_models, "No eligible models.", 10) })
  output$reference_table <- renderDT({ req(results()); safe_dt(results()$reference_summary, "No reference summary available.", 10) })
  output$third_drug_selection_table <- renderDT({ req(results()); safe_dt(results()$third_drug_selection, "No third drugs selected.", 25) })
  output$monotherapy_table <- renderDT({ req(results()); safe_dt(results()$monotherapy_summary, "No monotherapy response data available.", 25) })
  output$monotherapy_dose_table <- renderDT({ req(results()); safe_dt(results()$monotherapy_dose_response, "No monotherapy dose-response data available.", 25) })
  output$ranking_table <- renderDT({ req(results()); safe_dt(results()$third_ranking, "No third-drug ranking available.", 25) })
  output$top10_table <- renderDT({ req(results()); safe_dt(results()$top10, "No candidate drugs available.", 25) })
  output$functional_table <- renderDT({ req(results()); safe_dt(results()$functional_tbl, "No functional CRISPR dependency evidence available.", 25) })
  output$molecular_table <- renderDT({ req(results()); safe_dt(results()$molecular_tbl, "No molecular target-expression evidence available.", 25) })
  output$biomarker_table <- renderDT({ req(results()); safe_dt(results()$biomarker_tbl, "No optional mutation biomarker evidence available.", 25) })
  output$graph_ranking_table <- renderDT({
    req(results())
    safe_dt(results()$graph_ranking, "Graph learning was not enabled or no graph results were produced.", 25)
  })
  output$graph_diagnostics_table <- renderDT({
    req(results())
    safe_dt(results()$graph_diagnostics, "No graph diagnostics available.", 25)
  })
  output$validation_summary_table <- renderDT({
    req(results())
    safe_dt(results()$validation_summary, "Enable the graph model and Validation Agent, then rerun the analysis.", 25)
  })
  output$validation_fold_table <- renderDT({
    req(results())
    safe_dt(results()$validation_fold_metrics, "No cross-validation fold metrics are available.", 25)
  })
  output$clustering_validation_table <- renderDT({
    req(results())
    safe_dt(results()$validation_clustering, "No graph-clustering validation results are available.", 25)
  })

  output$monotherapy_plot <- renderPlot({
    req(results())
    plot_data <- results()$monotherapy_summary
    validate(need(nrow(plot_data) > 0, "No eligible monotherapy results are available for plotting."))

    plot_data <- plot_data %>%
      group_by(CellLineName) %>%
      slice_max(MeanMonotherapyResponse, n = 10, with_ties = FALSE) %>%
      ungroup()

    ggplot(plot_data, aes(x = reorder(ThirdDrug, MeanMonotherapyResponse), y = MeanMonotherapyResponse, fill = MonotherapyClass)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Top monotherapy responses in matched cervical models",
        subtitle = "Predicted inhibition is derived from the fitted GDSC2 response curve",
        x = "Third drug",
        y = "Mean predicted monotherapy inhibition (%)",
        fill = "Monotherapy class"
      )
  })

  output$monotherapy_dose_plot <- renderPlot({
    req(results())
    plot_data <- results()$monotherapy_dose_response
    validate(need(nrow(plot_data) > 0, "No monotherapy dose-response data are available for plotting."))

    top_ids <- results()$monotherapy_summary %>%
      group_by(CellLineName) %>%
      slice_max(MeanMonotherapyResponse, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      select(ModelID, ThirdDrugCompoundID)

    plot_data <- plot_data %>%
      semi_join(top_ids, by = c("ModelID", "ThirdDrugCompoundID"))

    ggplot(plot_data, aes(x = Dose, y = PredictedMonotherapyInhibition, group = ThirdDrug, linetype = ThirdDrug)) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.6) +
      scale_x_log10() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = "Third-drug monotherapy dose-response profiles",
        subtitle = "Only the top ten monotherapy candidates per model are shown",
        x = "Dose",
        y = "Predicted monotherapy inhibition (%)",
        linetype = "Third drug"
      )
  })

  output$top10_plot <- renderPlot({
    req(results())
    plot_data <- results()$top10
    validate(need(nrow(plot_data) > 0, "No eligible third-drug results are available for plotting."))

    reference_lines <- results()$reference_summary %>%
      distinct(ModelID, CellLineName, ReferenceMeanSynergy)

    ggplot(plot_data, aes(x = reorder(ThirdDrug, MeanSynergy), y = MeanSynergy, fill = Interaction)) +
      geom_col() +
      geom_hline(
        data = reference_lines,
        aes(yintercept = ReferenceMeanSynergy),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = paste0(
          "Cervical models: ",
          str_replace_all(input$combo, "_", " "),
          " + top third drugs (",
          synergy_method_label(input$synergy_method),
          ")"
        ),
        subtitle = paste(
          "Dashed line = model-specific reference mean;",
          "interaction model =", synergy_method_label(input$synergy_method)
        ),
        x = "Third drug",
        y = "Mean synergy score",
        fill = "Interaction"
      )
  })

  output$mono_combination_plot <- renderPlot({
    req(results())
    plot_data <- results()$top10
    validate(need(nrow(plot_data) > 0, "No monotherapy-versus-combination results are available."))

    ggplot(plot_data, aes(x = MeanMonotherapyResponse, y = MeanObserved, label = ThirdDrug, shape = Interaction)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_point(size = 2.4) +
      geom_text(check_overlap = TRUE, nudge_y = 2, size = 3) +
      facet_wrap(~ CellLineName, scales = "free") +
      theme_bw(base_size = 11) +
      labs(
        title = "Intrinsic monotherapy response versus observed combination response",
        subtitle = "Points above the diagonal show greater observed combination inhibition than monotherapy inhibition",
        x = "Mean predicted monotherapy inhibition (%)",
        y = "Mean observed combination inhibition (%)",
        shape = "Interaction"
      )
  })

  output$graph_ranking_plot <- renderPlot({
    req(results())
    plot_data <- results()$graph_ranking
    validate(need(nrow(plot_data) > 0, "Enable graph learning and rerun the analysis."))

    plot_data <- plot_data %>%
      group_by(CellLineName) %>%
      slice_max(GraphProbability, n = 10, with_ties = FALSE) %>%
      ungroup()

    ggplot(plot_data, aes(
      x = reorder(ThirdDrug, GraphProbability),
      y = GraphProbability,
      fill = GraphPredictedClass
    )) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Graph neural network prioritisation of candidate third drugs",
        subtitle = "Probability that the candidate exceeds the model-specific clinical reference",
        x = "Third drug",
        y = "Graph-predicted probability",
        fill = "Prediction"
      )
  })

  output$drug_graph_plot <- renderPlot({
    req(results(), input$graph_model_view)
    ranking <- results()$graph_ranking %>% filter(ModelID == input$graph_model_view)
    edges <- results()$graph_edges %>% filter(ModelID == input$graph_model_view)
    validate(need(nrow(ranking) > 0, "No graph data available for this model."))

    vertices <- ranking %>%
      transmute(
        name = ThirdDrug,
        probability = GraphProbability,
        label_class = GraphPredictedClass,
        observed_class = Interaction
      )

    if (nrow(edges) > 0) {
      g <- igraph::graph_from_data_frame(
        edges %>% select(from, to, weight),
        directed = FALSE,
        vertices = vertices
      )
      E(g)$width <- 0.5 + 4 * E(g)$weight
    } else {
      g <- igraph::make_empty_graph(n = nrow(vertices), directed = FALSE)
      V(g)$name <- vertices$name
      V(g)$probability <- vertices$probability
      V(g)$label_class <- vertices$label_class
    }

    vertex_size <- 6 + 14 * scales::rescale(V(g)$probability, to = c(0, 1))
    vertex_shape <- ifelse(V(g)$label_class == "Predicted stronger than reference", "circle", "square")

    plot(
      g,
      layout = igraph::layout_with_fr(g),
      vertex.size = vertex_size,
      vertex.shape = vertex_shape,
      vertex.label.cex = 0.65,
      vertex.label.dist = 0.4,
      edge.curved = 0.10,
      main = paste0(
        "Candidate-drug graph: ", ranking$CellLineName[[1]],
        " (", synergy_method_label(input$synergy_method), ")"
      )
    )
  })

  output$validation_roc_plot <- renderPlot({
    req(results())
    plot_data <- results()$validation_roc
    validate(need(nrow(plot_data) > 0, "No cross-validated ROC data are available."))

    ggplot(plot_data, aes(x = FalsePositiveRate, y = TruePositiveRate, linetype = Benchmark)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ CellLineName) +
      coord_equal() +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = "Cross-validated ROC curves for graph-model benchmarking",
        subtitle = paste("Synergy model:", synergy_method_label(input$synergy_method)),
        x = "False-positive rate", y = "True-positive rate", linetype = "Benchmark"
      )
  })

  output$validation_confusion_plot <- renderPlot({
    req(results())
    plot_data <- results()$validation_confusion
    validate(need(nrow(plot_data) > 0, "No cross-validated confusion-matrix data are available."))

    ggplot(plot_data, aes(x = Predicted, y = Observed, fill = Count)) +
      geom_tile() +
      geom_text(aes(label = Count), size = 4) +
      facet_wrap(~ Benchmark) +
      theme_bw(base_size = 11) +
      labs(
        title = paste0("Cross-validated confusion matrices at threshold ", input$validation_threshold),
        x = "Predicted class", y = "Observed class", fill = "Count"
      )
  })

  output$decision_curve_plot <- renderPlot({
    req(results())
    plot_data <- results()$validation_decision_curve
    validate(need(nrow(plot_data) > 0, "No decision-curve data are available."))

    reference <- plot_data %>%
      distinct(ModelID, CellLineName, Threshold, TreatAllNetBenefit, TreatNoneNetBenefit) %>%
      pivot_longer(
        cols = c(TreatAllNetBenefit, TreatNoneNetBenefit),
        names_to = "ReferenceStrategy", values_to = "NetBenefit"
      ) %>%
      mutate(ReferenceStrategy = recode(
        ReferenceStrategy,
        TreatAllNetBenefit = "Treat all",
        TreatNoneNetBenefit = "Treat none"
      ))

    ggplot() +
      geom_line(
        data = plot_data,
        aes(x = Threshold, y = NetBenefit, linetype = Benchmark),
        linewidth = 0.9
      ) +
      geom_line(
        data = reference,
        aes(x = Threshold, y = NetBenefit, linetype = ReferenceStrategy),
        linewidth = 0.7
      ) +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      labs(
        title = "Decision Curve Analysis of cross-validated graph predictions",
        subtitle = "Net benefit across clinically relevant probability thresholds",
        x = "Decision threshold", y = "Net benefit", linetype = "Model or strategy"
      )
  })

  output$clustering_validation_plot <- renderPlot({
    req(results())
    plot_data <- results()$validation_clustering %>%
      select(CellLineName, SilhouetteScore, RandIndex, AdjustedRandIndex) %>%
      pivot_longer(
        cols = c(SilhouetteScore, RandIndex, AdjustedRandIndex),
        names_to = "Metric", values_to = "Score"
      )
    validate(need(any(is.finite(plot_data$Score)), "No graph-clustering validation scores are available."))

    ggplot(plot_data, aes(x = CellLineName, y = Score, fill = Metric)) +
      geom_col(position = position_dodge(width = 0.8)) +
      coord_flip() +
      theme_bw(base_size = 11) +
      labs(
        title = "Graph-clustering validity and label agreement",
        subtitle = "Silhouette score measures graph-cluster separation; Rand indices measure agreement with interaction labels",
        x = "Cervical cancer model", y = "Validation score", fill = "Metric"
      )
  })

  output$functional_plot <- renderPlot({
    req(results())
    plot_data <- results()$functional_tbl %>%
      filter(!is.na(CRISPR_Dependency), DependencyClass %in% c("Strong dependency", "Moderate dependency", "Weak dependency"))

    validate(need(nrow(plot_data) > 0, "No CRISPR dependency data available for plotted targets."))

    ggplot(plot_data, aes(x = reorder(TargetGene, CRISPR_Dependency), y = CRISPR_Dependency, fill = DependencyClass)) +
      geom_col() +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_hline(yintercept = 0.8, linetype = "dotted") +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Functional evidence: CRISPR dependency of third-drug targets",
        subtitle = "Dashed = moderate dependency cutoff; dotted = strong dependency cutoff",
        x = "Target gene",
        y = "CRISPR dependency score",
        fill = "Dependency class"
      )
  })

  output$molecular_plot <- renderPlot({
    req(results())
    plot_data <- results()$molecular_tbl %>% filter(!is.na(Expression_log2TPM))
    validate(need(nrow(plot_data) > 0, "No expression data available for plotted targets."))

    ggplot(plot_data, aes(x = reorder(TargetGene, Expression_log2TPM), y = Expression_log2TPM, fill = ExpressionClass)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ CellLineName, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(
        title = "Molecular evidence: target expression in matched cervical models",
        subtitle = "Expression values are log2(TPM + 1); expression classes use target-specific percentiles",
        x = "Target gene",
        y = "log2(TPM + 1)",
        fill = "Expression class"
      )
  })

  output$biomarker_plot <- renderPlot({
    req(results())
    plot_data <- results()$biomarker_tbl %>% count(CellLineName, MutationStatus, name = "n")
    validate(need(nrow(plot_data) > 0, "No mutation biomarker data available."))

    ggplot(plot_data, aes(x = CellLineName, y = n, fill = MutationStatus)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 11) +
      labs(
        title = "Optional biomarker evidence: target mutation context",
        x = "Cervical cancer model",
        y = "Drug-target associations",
        fill = "Mutation status"
      )
  })

  output$download_antiangiogenic_set <- downloadHandler(
    filename = function() {
      "Cervical_GDSC2_Antiangiogenic_Compound_Set.csv"
    },
    content = function(file) {
      write.csv(
        Cervical_Antiangiogenic_Compounds,
        file,
        row.names = FALSE
      )
    }
  )

  output$download_reference <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_Reference_Summary.csv"
    ),
    content = function(file) write.csv(results()$reference_summary, file, row.names = FALSE)
  )

  output$download_monotherapy <- downloadHandler(
    filename = function() paste0(input$combo, "_Monotherapy_Response_Summary.csv"),
    content = function(file) write.csv(results()$monotherapy_summary, file, row.names = FALSE)
  )

  output$download_monotherapy_dose <- downloadHandler(
    filename = function() paste0(input$combo, "_Monotherapy_Dose_Response.csv"),
    content = function(file) write.csv(results()$monotherapy_dose_response, file, row.names = FALSE)
  )

  output$download_ranking <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_ThirdDrug_Ranking.csv"
    ),
    content = function(file) write.csv(results()$third_ranking, file, row.names = FALSE)
  )

  output$download_top10 <- downloadHandler(
    filename = function() paste0(
      input$combo, "_", toupper(input$synergy_method), "_Top_Candidates_Per_Model.csv"
    ),
    content = function(file) write.csv(results()$top10, file, row.names = FALSE)
  )

  output$download_functional <- downloadHandler(
    filename = function() paste0(input$combo, "_Functional_CRISPR_Evidence.csv"),
    content = function(file) write.csv(results()$functional_tbl, file, row.names = FALSE)
  )

  output$download_molecular <- downloadHandler(
    filename = function() paste0(input$combo, "_Molecular_Expression_Evidence.csv"),
    content = function(file) write.csv(results()$molecular_tbl, file, row.names = FALSE)
  )

  output$download_biomarker <- downloadHandler(
    filename = function() paste0(input$combo, "_Optional_Mutation_Biomarker_Evidence.csv"),
    content = function(file) write.csv(results()$biomarker_tbl, file, row.names = FALSE)
  )

  output$download_graph_ranking <- downloadHandler(
    filename = function() paste0(input$combo, "_", toupper(input$synergy_method), "_GraphML_Ranking.csv"),
    content = function(file) write.csv(results()$graph_ranking, file, row.names = FALSE)
  )

  output$download_graph_edges <- downloadHandler(
    filename = function() paste0(input$combo, "_Drug_Graph_Edges.csv"),
    content = function(file) write.csv(results()$graph_edges, file, row.names = FALSE)
  )

  output$download_graph_diagnostics <- downloadHandler(
    filename = function() paste0(input$combo, "_GraphML_Diagnostics.csv"),
    content = function(file) write.csv(results()$graph_diagnostics, file, row.names = FALSE)
  )

  output$download_validation_summary <- downloadHandler(
    filename = function() paste0(input$combo, "_", toupper(input$synergy_method), "_Validation_Agent_Summary.csv"),
    content = function(file) write.csv(results()$validation_summary, file, row.names = FALSE)
  )

  output$download_validation_predictions <- downloadHandler(
    filename = function() paste0(input$combo, "_", toupper(input$synergy_method), "_CrossValidation_Predictions.csv"),
    content = function(file) write.csv(results()$validation_predictions, file, row.names = FALSE)
  )

  output$download_validation_folds <- downloadHandler(
    filename = function() paste0(input$combo, "_", toupper(input$synergy_method), "_CrossValidation_Fold_Metrics.csv"),
    content = function(file) write.csv(results()$validation_fold_metrics, file, row.names = FALSE)
  )

  output$download_decision_curve <- downloadHandler(
    filename = function() paste0(input$combo, "_", toupper(input$synergy_method), "_Decision_Curve_Data.csv"),
    content = function(file) write.csv(results()$validation_decision_curve, file, row.names = FALSE)
  )
}
