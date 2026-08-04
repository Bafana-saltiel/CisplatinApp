# ============================================================
# 10. CROSS-SCC COHORT VALIDATION MODULE
# ============================================================

CROSS_SCC_DATA_DIR <- file.path("data", "TCGA data")

CROSS_SCC_TARGETS <- c(
  "CDK2", "DOT1L", "TOP1", "TUBA1B", "BRD1", "PBRM1",
  "PARP2", "ATR", "EWSR1", "BRAF", "MTOR"
)

CROSS_SCC_PROGRAMMES <- list(
  "Multi-target DDR-epigenetic-mTOR" = c("ATR", "PARP2", "DOT1L", "MTOR"),
  "DNA damage and replication" = c("ATR", "PARP2", "TOP1"),
  "Chromatin and transcription" = c("DOT1L", "BRD1", "PBRM1", "EWSR1"),
  "Cell cycle" = "CDK2",
  "Microtubule organisation" = "TUBA1B",
  "MAPK signalling" = "BRAF",
  "mTOR signalling" = "MTOR"
)

cross_scc_required_files <- function() {
  projects <- c("CESC", "HNSC", "LUSC", "ESCA")
  c(
    file.path(CROSS_SCC_DATA_DIR, "GSE299125_expression_matrix.csv"),
    file.path(CROSS_SCC_DATA_DIR, "GSE299125_sample_metadata.csv"),
    unlist(lapply(projects, function(project) {
      file.path(
        CROSS_SCC_DATA_DIR,
        c(
          paste0("TCGA_", project, "_expression_matrix.csv"),
          paste0("TCGA_", project, "_sample_metadata.csv")
        )
      )
    }))
  )
}

cross_scc_safe_z <- function(x) {
  x <- as.numeric(x)
  value_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / value_sd
}

cross_scc_is_scc <- function(primary_diagnosis, morphology) {
  diagnosis <- stringr::str_to_lower(dplyr::coalesce(primary_diagnosis, ""))
  diagnosis_flag <- stringr::str_detect(diagnosis, "squamous|basaloid") &
    !stringr::str_detect(diagnosis, "adenosquamous|adenocarcinoma")

  morphology_number <- suppressWarnings(as.integer(
    stringr::str_extract(as.character(morphology), "^[0-9]{4}")
  ))
  morphology_flag <- !is.na(morphology_number) &
    morphology_number >= 8070 & morphology_number <= 8084

  diagnosis_flag | morphology_flag
}

cross_scc_validate_expression <- function(expression_data, label) {
  required_columns <- c("GeneSymbol", "EnsemblGeneID")
  missing_columns <- setdiff(required_columns, names(expression_data))
  if (length(missing_columns) > 0) {
    stop(label, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }

  missing_targets <- setdiff(CROSS_SCC_TARGETS, expression_data$GeneSymbol)
  if (length(missing_targets) > 0) {
    stop(label, " is missing targets: ", paste(missing_targets, collapse = ", "))
  }

  if (anyDuplicated(expression_data$GeneSymbol)) {
    stop(label, " contains duplicated gene symbols.")
  }

  expression_data %>%
    filter(GeneSymbol %in% CROSS_SCC_TARGETS) %>%
    arrange(match(GeneSymbol, CROSS_SCC_TARGETS))
}

cross_scc_make_matrix <- function(expression_data, sample_ids, label) {
  absent_samples <- setdiff(sample_ids, names(expression_data))
  if (length(absent_samples) > 0) {
    stop(label, " is missing samples: ", paste(head(absent_samples, 10), collapse = ", "))
  }

  matrix_data <- expression_data %>%
    select(GeneSymbol, all_of(sample_ids)) %>%
    tibble::column_to_rownames("GeneSymbol") %>%
    as.matrix()
  storage.mode(matrix_data) <- "numeric"

  if (any(!is.finite(matrix_data))) {
    stop(label, " contains missing or non-finite expression values.")
  }
  matrix_data
}

cross_scc_load_gse <- function() {
  expression_data <- readr::read_csv(
    file.path(CROSS_SCC_DATA_DIR, "GSE299125_expression_matrix.csv"),
    show_col_types = FALSE
  ) %>% cross_scc_validate_expression("GSE299125 expression")

  metadata <- readr::read_csv(
    file.path(CROSS_SCC_DATA_DIR, "GSE299125_sample_metadata.csv"),
    show_col_types = FALSE
  )

  required_metadata <- c("PatientID", "HIV_Status", "Disease_Status")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0) {
    stop("GSE299125 metadata is missing: ", paste(missing_metadata, collapse = ", "))
  }

  sample_ids <- setdiff(names(expression_data), c("GeneSymbol", "EnsemblGeneID"))
  metadata <- metadata %>%
    filter(PatientID %in% sample_ids) %>%
    slice(match(sample_ids, PatientID)) %>%
    mutate(
      Cohort = "GSE299125",
      SCC = Disease_Status == "Invasive cervical carcinoma",
      IncludePrimaryAnalysis = SCC,
      HistologyGroup = Disease_Status
    )

  if (any(is.na(metadata$PatientID))) {
    stop("Not all GSE299125 expression columns could be matched to metadata.")
  }

  raw_matrix <- cross_scc_make_matrix(expression_data, sample_ids, "GSE299125")

  # Whole-transcriptome library sizes are unavailable in the reduced 11-gene
  # file. log2(count + 1) is used only to stabilize variance before calculating
  # within-gene, within-cohort Z-scores. Absolute GSE and TCGA values are never
  # compared directly.
  transformed <- log2(raw_matrix + 1)
  invasive_ids <- metadata$PatientID[metadata$IncludePrimaryAnalysis]
  invasive_matrix <- transformed[, invasive_ids, drop = FALSE]
  z_matrix <- t(apply(invasive_matrix, 1, cross_scc_safe_z))
  rownames(z_matrix) <- rownames(invasive_matrix)
  colnames(z_matrix) <- colnames(invasive_matrix)

  list(
    cohort = "GSE299125",
    expression = invasive_matrix,
    z = z_matrix,
    metadata_all = metadata,
    metadata = metadata %>% filter(IncludePrimaryAnalysis)
  )
}

cross_scc_load_tcga <- function(project) {
  cohort <- paste0("TCGA-", project)
  expression_data <- readr::read_csv(
    file.path(CROSS_SCC_DATA_DIR, paste0("TCGA_", project, "_expression_matrix.csv")),
    show_col_types = FALSE
  ) %>% cross_scc_validate_expression(paste0(cohort, " expression"))

  metadata <- readr::read_csv(
    file.path(CROSS_SCC_DATA_DIR, paste0("TCGA_", project, "_sample_metadata.csv")),
    show_col_types = FALSE
  )

  required_metadata <- c(
    "PatientID", "SampleType", "PrimaryDiagnosis", "Morphology",
    "TumorStage", "OverallSurvivalDays", "OverallSurvivalStatus"
  )
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0) {
    stop(cohort, " metadata is missing: ", paste(missing_metadata, collapse = ", "))
  }

  sample_ids <- setdiff(names(expression_data), c("GeneSymbol", "EnsemblGeneID"))
  metadata <- metadata %>%
    filter(PatientID %in% sample_ids) %>%
    slice(match(sample_ids, PatientID)) %>%
    mutate(
      Cohort = cohort,
      PrimaryTumour = str_to_lower(SampleType) == "primary tumor",
      SCC = cross_scc_is_scc(PrimaryDiagnosis, Morphology),
      IncludePrimaryAnalysis = PrimaryTumour & SCC,
      HistologyGroup = if_else(SCC, "Squamous cell carcinoma", "Non-SCC")
    )

  if (any(is.na(metadata$PatientID))) {
    stop("Not all ", cohort, " expression columns could be matched to metadata.")
  }

  included_ids <- metadata$PatientID[metadata$IncludePrimaryAnalysis]
  if (length(included_ids) < 2) {
    stop(cohort, " has fewer than two eligible primary SCC samples.")
  }

  expression_matrix <- cross_scc_make_matrix(expression_data, included_ids, cohort)
  z_matrix <- t(apply(expression_matrix, 1, cross_scc_safe_z))
  rownames(z_matrix) <- rownames(expression_matrix)
  colnames(z_matrix) <- colnames(expression_matrix)

  list(
    cohort = cohort,
    expression = expression_matrix,
    z = z_matrix,
    metadata_all = metadata,
    metadata = metadata %>% filter(IncludePrimaryAnalysis)
  )
}

cross_scc_load_all <- function() {
  missing_files <- cross_scc_required_files()[!file.exists(cross_scc_required_files())]
  if (length(missing_files) > 0) {
    stop("Missing cross-SCC files: ", paste(missing_files, collapse = ", "))
  }

  cohorts <- c(
    list(GSE299125 = cross_scc_load_gse()),
    setNames(
      lapply(c("CESC", "HNSC", "LUSC", "ESCA"), cross_scc_load_tcga),
      paste0("TCGA-", c("CESC", "HNSC", "LUSC", "ESCA"))
    )
  )

  z_long <- purrr::imap_dfr(cohorts, function(cohort_data, cohort_name) {
    as.data.frame(cohort_data$z, check.names = FALSE) %>%
      rownames_to_column("GeneSymbol") %>%
      pivot_longer(-GeneSymbol, names_to = "PatientID", values_to = "ZScore") %>%
      mutate(Cohort = cohort_name, .before = 1)
  })

  metadata <- bind_rows(lapply(cohorts, function(x) x$metadata))

  audit <- bind_rows(lapply(cohorts, function(x) {
    x$metadata_all %>%
      count(Cohort, HistologyGroup, IncludePrimaryAnalysis, name = "N")
  })) %>% arrange(Cohort, desc(IncludePrimaryAnalysis), HistologyGroup)

  list(cohorts = cohorts, z_long = z_long, metadata = metadata, audit = audit)
}

cross_scc_programme_scores <- function(analysis_data) {
  purrr::imap_dfr(analysis_data$cohorts, function(cohort_data, cohort_name) {
    purrr::map_dfr(names(CROSS_SCC_PROGRAMMES), function(programme_name) {
      genes <- CROSS_SCC_PROGRAMMES[[programme_name]]
      tibble(
        Cohort = cohort_name,
        PatientID = colnames(cohort_data$z),
        Programme = programme_name,
        Score = colMeans(cohort_data$z[genes, , drop = FALSE])
      )
    })
  })
}

cross_scc_state_projection <- function(analysis_data, k) {
  gse_z <- analysis_data$cohorts$GSE299125$z
  patient_matrix <- t(gse_z)
  tree <- hclust(dist(patient_matrix), method = "ward.D2")
  discovery_state <- cutree(tree, k = k)

  centroids <- sapply(sort(unique(discovery_state)), function(state) {
    rowMeans(gse_z[, discovery_state == state, drop = FALSE])
  })
  colnames(centroids) <- paste0("State ", sort(unique(discovery_state)))

  discovery <- tibble(
    Cohort = "GSE299125",
    PatientID = names(discovery_state),
    AssignedState = paste0("State ", unname(discovery_state)),
    DistanceToCentroid = NA_real_,
    AssignmentMargin = NA_real_,
    Projection = "Discovery"
  )

  projected <- purrr::imap_dfr(
    analysis_data$cohorts[names(analysis_data$cohorts) != "GSE299125"],
    function(cohort_data, cohort_name) {
      samples <- t(cohort_data$z[rownames(centroids), , drop = FALSE])
      distances <- sapply(seq_len(ncol(centroids)), function(index) {
        centroid_matrix <- matrix(
          centroids[, index], nrow = nrow(samples),
          ncol = ncol(samples), byrow = TRUE
        )
        sqrt(rowSums((samples - centroid_matrix)^2))
      })
      colnames(distances) <- colnames(centroids)
      nearest <- max.col(-distances, ties.method = "first")
      sorted_distances <- t(apply(distances, 1, sort))

      tibble(
        Cohort = cohort_name,
        PatientID = rownames(samples),
        AssignedState = colnames(centroids)[nearest],
        DistanceToCentroid = distances[cbind(seq_len(nrow(distances)), nearest)],
        AssignmentMargin = sorted_distances[, 2] - sorted_distances[, 1],
        Projection = "External projection"
      )
    }
  )

  list(assignments = bind_rows(discovery, projected), centroids = centroids)
}

cross_scc_network_preservation <- function(analysis_data) {
  upper_values <- function(correlation_matrix) {
    correlation_matrix[upper.tri(correlation_matrix, diag = FALSE)]
  }

  correlation_matrices <- lapply(analysis_data$cohorts, function(x) {
    cor(t(x$z), method = "spearman", use = "pairwise.complete.obs")
  })
  reference_edges <- upper_values(correlation_matrices$GSE299125)

  purrr::imap_dfr(
    correlation_matrices[names(correlation_matrices) != "GSE299125"],
    function(validation_matrix, cohort_name) {
      validation_edges <- upper_values(validation_matrix)
      test <- suppressWarnings(cor.test(
        reference_edges, validation_edges,
        method = "spearman", exact = FALSE
      ))
      tibble(
        Reference = "GSE299125 invasive cervical carcinoma",
        ValidationCohort = cohort_name,
        TargetPairs = length(reference_edges),
        NetworkConcordance = unname(test$estimate),
        PValue = test$p.value
      )
    }
  ) %>% mutate(FDR = p.adjust(PValue, method = "BH"))
}

cross_scc_validation_ui <- function() {
  tagList(
    fluidRow(
      column(
        width = 4,
        wellPanel(
          h4("Analysis controls"),
          selectInput(
            "cross_scc_cohort",
            "Cohort displayed in the patient heatmap:",
            choices = c("GSE299125", "TCGA-CESC", "TCGA-HNSC", "TCGA-LUSC", "TCGA-ESCA"),
            selected = "GSE299125"
          ),
          sliderInput(
            "cross_scc_high_threshold",
            "Target-high Z-score threshold:",
            min = 0.5, max = 2.0, value = 1.0, step = 0.1
          ),
          selectInput(
            "cross_scc_state_k",
            "Number of GSE299125-derived states:",
            choices = c("2 states" = 2, "3 states" = 3, "4 states" = 4),
            selected = 3
          ),
          actionButton(
            "run_cross_scc",
            "Run cross-SCC validation",
            class = "btn-success",
            icon = icon("play")
          ),
          br(), br(),
          downloadButton("download_cross_scc_table", "Download displayed table")
        )
      ),
      column(
        width = 8,
        wellPanel(
          h4("Interpretation"),
          p(
            "The fixed 11-target panel is derived from cervical model-specific pharmacology and CRISPR dependency. GSE299125 provides independent cervical patient-level molecular validation; TCGA-CESC tests within-disease replication; and TCGA-HNSC, TCGA-LUSC and TCGA-ESCA-SCC test cross-squamous conservation."
          ),
          p(
            "GSE299125 raw counts and TCGA log2(FPKM-UQ + 1) values are not compared directly. Each target is standardized within its cohort, so results represent relative patient states rather than absolute cross-platform abundance or predicted clinical response."
          ),
          uiOutput("cross_scc_status")
        )
      )
    ),
    tabsetPanel(
      id = "cross_scc_results_tab",
      tabPanel(
        "Cohort inclusion",
        h4("Auditable histology and primary-tumour filtering"),
        DTOutput("cross_scc_audit_table")
      ),
      tabPanel(
        "Patient target states",
        plotOutput("cross_scc_heatmap", height = "800px")
      ),
      tabPanel(
        "Target-high prevalence",
        plotOutput("cross_scc_prevalence_plot", height = "700px"),
        DTOutput("cross_scc_prevalence_table")
      ),
      tabPanel(
        "Vulnerability programmes",
        plotOutput("cross_scc_programme_plot", height = "800px")
      ),
      tabPanel(
        "Projected states",
        plotOutput("cross_scc_state_plot", height = "650px"),
        DTOutput("cross_scc_state_table")
      ),
      tabPanel(
        "Cross-SCC PCA",
        plotOutput("cross_scc_pca_plot", height = "700px")
      ),
      tabPanel(
        "Network preservation",
        h4("Preservation of the 55 pairwise target correlations"),
        p("Interpret cautiously because the GSE299125 invasive reference contains 13 patients."),
        DTOutput("cross_scc_network_table")
      )
    )
  )
}

cross_scc_validation_server <- function(input, output, session) {
  analysis_data <- eventReactive(input$run_cross_scc, {
    withProgress(message = "Loading and harmonizing validation cohorts...", value = 0, {
      incProgress(0.2, detail = "Checking expression and metadata files")
      data <- cross_scc_load_all()
      incProgress(0.8, detail = "Preparing cross-cohort target states")
      data
    })
  }, ignoreInit = TRUE)

  output$cross_scc_status <- renderUI({
    missing_files <- cross_scc_required_files()[!file.exists(cross_scc_required_files())]
    if (length(missing_files) > 0) {
      return(tags$div(
        class = "alert alert-danger",
        tags$b("Cross-SCC validation files are incomplete."),
        tags$br(),
        paste(basename(missing_files), collapse = ", ")
      ))
    }

    if (input$run_cross_scc == 0) {
      return(tags$div(
        class = "alert alert-info",
        "All 10 expression/metadata files were detected. Click Run cross-SCC validation."
      ))
    }

    data <- analysis_data()
    counts <- vapply(data$cohorts, function(x) ncol(x$z), integer(1))
    tags$div(
      class = "alert alert-success",
      tags$b("Validation cohorts loaded: "),
      paste(paste(names(counts), counts, sep = " = "), collapse = "; ")
    )
  })

  output$cross_scc_audit_table <- renderDT({
    req(analysis_data())
    datatable(
      analysis_data()$audit,
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  target_prevalence <- reactive({
    req(analysis_data())
    analysis_data()$z_long %>%
      group_by(Cohort, GeneSymbol) %>%
      summarise(
        N = n(),
        TargetHighN = sum(ZScore >= input$cross_scc_high_threshold),
        TargetHighProportion = mean(ZScore >= input$cross_scc_high_threshold),
        MedianZScore = median(ZScore),
        .groups = "drop"
      )
  })

  output$cross_scc_prevalence_table <- renderDT({
    datatable(
      target_prevalence() %>%
        mutate(TargetHighProportion = round(TargetHighProportion, 3)),
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  output$cross_scc_prevalence_plot <- renderPlot({
    target_prevalence() %>%
      mutate(GeneSymbol = factor(GeneSymbol, levels = CROSS_SCC_TARGETS)) %>%
      ggplot(aes(GeneSymbol, TargetHighProportion, fill = Cohort)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.75) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = NULL,
        y = paste0("Proportion with Z-score >= ", input$cross_scc_high_threshold),
        fill = "Cohort",
        title = "Prevalence of comparatively elevated therapeutic-vulnerability targets"
      ) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$cross_scc_heatmap <- renderPlot({
    req(analysis_data(), input$cross_scc_cohort)
    cohort_data <- analysis_data()$cohorts[[input$cross_scc_cohort]]
    validate(need(!is.null(cohort_data), "The selected cohort is unavailable."))

    heatmap_data <- as.data.frame(cohort_data$z, check.names = FALSE) %>%
      rownames_to_column("GeneSymbol") %>%
      pivot_longer(-GeneSymbol, names_to = "PatientID", values_to = "ZScore") %>%
      mutate(
        GeneSymbol = factor(GeneSymbol, levels = rev(CROSS_SCC_TARGETS)),
        PatientID = factor(PatientID, levels = colnames(cohort_data$z))
      )

    ggplot(heatmap_data, aes(PatientID, GeneSymbol, fill = pmax(-3, pmin(3, ZScore)))) +
      geom_tile() +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-3, 3), name = "Z-score"
      ) +
      labs(
        x = "Patient",
        y = NULL,
        title = paste0(input$cross_scc_cohort, ": relative patient target-expression states")
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 90, hjust = 1, vjust = 0.5,
          size = if (ncol(cohort_data$z) > 150) 3 else 6
        )
      )
  })

  programme_scores <- reactive({
    req(analysis_data())
    cross_scc_programme_scores(analysis_data())
  })

  output$cross_scc_programme_plot <- renderPlot({
    ggplot(programme_scores(), aes(Cohort, Score, fill = Cohort)) +
      geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.8) +
      geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
      facet_wrap(~Programme, ncol = 3, scales = "free_y") +
      labs(
        x = NULL,
        y = "Mean within-cohort target Z-score",
        title = "Model-derived vulnerability-programme scores"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })

  state_projection <- reactive({
    req(analysis_data())
    cross_scc_state_projection(
      analysis_data(),
      k = as.integer(input$cross_scc_state_k)
    )
  })

  state_prevalence <- reactive({
    state_projection()$assignments %>%
      count(Cohort, AssignedState, name = "N") %>%
      group_by(Cohort) %>%
      mutate(Proportion = N / sum(N)) %>%
      ungroup()
  })

  output$cross_scc_state_plot <- renderPlot({
    ggplot(state_prevalence(), aes(Cohort, Proportion, fill = AssignedState)) +
      geom_col(width = 0.75) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = NULL, y = "Proportion of tumours", fill = "Projected state",
        title = "Projection of frozen GSE299125-derived states into TCGA SCC cohorts"
      ) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  output$cross_scc_state_table <- renderDT({
    datatable(
      state_projection()$assignments,
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  pca_results <- reactive({
    req(analysis_data())
    combined <- bind_rows(purrr::imap(analysis_data()$cohorts, function(x, cohort_name) {
      as.data.frame(t(x$z), check.names = FALSE) %>%
        rownames_to_column("PatientID") %>%
        mutate(Cohort = cohort_name, .after = PatientID)
    }))
    fit <- prcomp(combined[, CROSS_SCC_TARGETS], center = FALSE, scale. = FALSE)
    variance <- 100 * summary(fit)$importance[2, 1:2]
    scores <- bind_cols(
      combined[, c("PatientID", "Cohort")],
      as.data.frame(fit$x[, 1:2, drop = FALSE])
    )
    list(scores = scores, variance = variance)
  })

  output$cross_scc_pca_plot <- renderPlot({
    pca <- pca_results()
    ggplot(pca$scores, aes(PC1, PC2, colour = Cohort)) +
      geom_point(alpha = 0.65, size = 1.8) +
      labs(
        x = sprintf("PC1 (%.1f%%)", pca$variance[1]),
        y = sprintf("PC2 (%.1f%%)", pca$variance[2]),
        title = "Cross-SCC relative vulnerability-state space"
      ) +
      theme_bw(base_size = 12)
  })

  network_results <- reactive({
    req(analysis_data())
    cross_scc_network_preservation(analysis_data())
  })

  output$cross_scc_network_table <- renderDT({
    datatable(
      network_results() %>%
        mutate(across(c(NetworkConcordance, PValue, FDR), ~signif(.x, 4))),
      rownames = FALSE,
      options = list(dom = "t", scrollX = TRUE)
    )
  })

  displayed_download <- reactive({
    req(analysis_data())
    selected_tab <- input$cross_scc_results_tab %||% "Cohort inclusion"
    switch(
      selected_tab,
      "Cohort inclusion" = analysis_data()$audit,
      "Patient target states" = analysis_data()$z_long %>%
        filter(Cohort == input$cross_scc_cohort),
      "Target-high prevalence" = target_prevalence(),
      "Vulnerability programmes" = programme_scores(),
      "Projected states" = state_projection()$assignments,
      "Cross-SCC PCA" = pca_results()$scores,
      "Network preservation" = network_results(),
      analysis_data()$audit
    )
  })

  output$download_cross_scc_table <- downloadHandler(
    filename = function() {
      tab_name <- input$cross_scc_results_tab %||% "Cohort inclusion"
      safe_name <- str_replace_all(str_to_lower(tab_name), "[^a-z0-9]+", "_")
      paste0("Cross_SCC_", safe_name, ".csv")
    },
    content = function(file) {
      readr::write_csv(displayed_download(), file, na = "")
    }
  )
}

