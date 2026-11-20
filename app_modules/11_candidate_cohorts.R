candidate_cohort_ui <- function() {
  tagList(
    fluidRow(
      column(4, selectInput("candidate_cohort", "Cohort",
        c("GSE299125", "TCGA-CESC", "TCGA-ESCA", "TCGA-HNSC", "TCGA-LUSC"))),
      column(4, selectInput("candidate_cohort_groups", "Sample groups",
        c("Invasive carcinoma by HIV status" = "invasive",
          "All cervical samples by HIV and disease status" = "all"))),
      column(4, numericInput("candidate_cohort_top", "Include ranks 1 through N per cell line",
        value = 1, min = 1, max = 10))
    ),
    checkboxInput("candidate_cohort_strongest", "Show one target per candidate and model",
                  value = FALSE),
    p("Run Calculate expectations after changing the regimen, candidates or doses. Rows use the candidate ranks, annotated targets and CRISPR evidence from Target context for the completed screen. Tumour expression is descriptive."),
    textOutput("candidate_cohort_coverage"),
    div(style = "overflow-x: auto;",
        div(style = "min-width: 1150px;", plotOutput("candidate_cohort_heatmap", height = "auto"))),
    downloadButton("candidate_cohort_download", "Download expression summary"),
    DTOutput("candidate_cohort_summary")
  )
}

candidate_cohort_data <- function(cohort, grouping) {
  prefix <- gsub("-", "_", cohort, fixed = TRUE)
  expression <- read_csv(
    file.path(CROSS_SCC_DATA_DIR, paste0(prefix, "_candidate_expression.csv")),
    show_col_types = FALSE
  )
  metadata <- read_csv(
    file.path(CROSS_SCC_DATA_DIR, paste0(prefix, "_sample_metadata.csv")),
    show_col_types = FALSE
  ) %>% filter(PatientID %in% names(expression))

  if (cohort == "GSE299125") {
    if (grouping == "invasive") {
      metadata <- filter(metadata, Disease_Status == "Invasive cervical carcinoma")
    }
    metadata <- metadata %>% mutate(Group = if (grouping == "invasive") HIV_Status else
      Integrated_Patient_Group)
    unit <- "CPM"
  } else {
    metadata <- metadata %>% filter(SampleType == "Primary Tumor")
    if (grouping == "histology") {
      metadata <- metadata %>%
        mutate(Group = if_else(cross_scc_is_scc(PrimaryDiagnosis, Morphology), "SCC", "Non-SCC"))
    } else {
      metadata <- metadata %>% mutate(Group = case_when(
        str_detect(TumorStage, "^Stage (III|IV)") ~ "Stage III–IV",
        str_detect(TumorStage, "^Stage (II|I)") ~ "Stage I–II"
      ))
    }
    unit <- "log2(FPKM-UQ + 1)"
  }

  metadata <- metadata %>% filter(!is.na(Group)) %>% arrange(Group, PatientID)
  values <- cross_scc_make_matrix(expression, metadata$PatientID, cohort)
  list(values = values, metadata = metadata, unit = unit, cohort = cohort)
}

candidate_cohort_projection <- function(target_context, cohort, strongest) {
  annotated <- target_context %>%
    distinct(ModelID, ThirdDrugCompoundID, TargetGene, .keep_all = TRUE)
  rows <- annotated %>% filter(TargetGene %in% rownames(cohort$values))
  omitted <- nrow(annotated) - nrow(rows)
  validate(need(nrow(rows) > 0L, paste(
    "The displayed candidate ranks have no annotated targets with expression data in this cohort.",
    "Increase the number of ranks to check the next-ranked candidates."
  )))

  if (strongest) {
    rows <- rows %>% arrange(desc(CRISPR_Dependency), TargetGene) %>%
      group_by(ModelID, ThirdDrugCompoundID) %>% slice_head(n = 1) %>% ungroup()
  }
  rows <- rows %>% mutate(
    Context = paste(coalesce(CellLineName, ModelID),
                    coalesce(CompoundName, ThirdDrugCompoundID), TargetGene, sep = " | ")
  )
  values <- cohort$values[unique(rows$TargetGene), , drop = FALSE]
  transformed <- if (cohort$unit == "CPM") log2(values + 1) else values
  z <- transformed
  for (i in seq_len(nrow(z))) z[i, ] <- cross_scc_safe_z(transformed[i, ])
  names <- cohort$metadata$PatientID
  groups <- cohort$metadata$Group
  order_samples <- unlist(lapply(split(names, groups), function(samples) {
    if (length(samples) < 3L) return(samples)
    samples[hclust(dist(t(z[, samples, drop = FALSE])))$order]
  }), use.names = FALSE)

  summary <- as.data.frame(values) %>% rownames_to_column("TargetGene") %>%
    pivot_longer(-TargetGene, names_to = "PatientID", values_to = "Expression") %>%
    left_join(select(cohort$metadata, PatientID, Group), by = "PatientID") %>%
    group_by(TargetGene, Group) %>%
    summarise(Samples = n(), Detected = sum(Expression > 0),
              MedianExpression = median(Expression), .groups = "drop")

  list(rows = rows, z = z, omitted = omitted, unit = cohort$unit,
       cohort = cohort$cohort, metadata = cohort$metadata, samples = order_samples,
       summary = inner_join(rows, summary, by = "TargetGene", relationship = "many-to-many") %>%
         mutate(Cohort = cohort$cohort, ExpressionUnit = cohort$unit))
}

candidate_cohort_server <- function(input, output, session, ranked_target_context) {
  observeEvent(input$candidate_cohort, {
    choices <- if (input$candidate_cohort == "GSE299125") {
      c("Invasive carcinoma by HIV status" = "invasive",
        "All cervical samples by HIV and disease status" = "all")
    } else {
      c("Stages I–II and III–IV" = "stage", "Primary tumour histology" = "histology")
    }
    updateSelectInput(session, "candidate_cohort_groups", choices = choices)
  })
  cohort <- reactive({
    req(input$candidate_cohort, input$candidate_cohort_groups)
    grouping <- input$candidate_cohort_groups
    if (input$candidate_cohort == "GSE299125") req(grouping %in% c("invasive", "all"))
    else req(grouping %in% c("stage", "histology"))
    candidate_cohort_data(input$candidate_cohort, grouping)
  })
  projection <- reactive({
    req(input$candidate_cohort_top >= 1, input$candidate_cohort_top <= 10)
    targets <- ranked_target_context() %>%
      filter(CandidateRank <= input$candidate_cohort_top)
    candidate_cohort_projection(targets, cohort(), isTRUE(input$candidate_cohort_strongest))
  })
  output$candidate_cohort_coverage <- renderText({
    p <- projection()
    paste(nrow(p$rows), "model–compound–target rows,",
          n_distinct(p$rows$TargetGene), "unique genes and", length(p$samples), "samples.",
          p$omitted, "annotations excluded for absent gene symbols or cohort expression.")
  })
  output$candidate_cohort_heatmap <- renderPlot({
    p <- projection()
    rows <- p$rows
    z <- p$z[rows$TargetGene, p$samples, drop = FALSE]
    order <- if (nrow(z) > 1) hclust(dist(z))$order else 1L
    rownames(z) <- rows$Context
    labels <- setNames(paste(rows$Context, rows$DependencyClass, sep = "  ·  "), rows$Context)
    heat <- as.data.frame(z, check.names = FALSE) %>%
      rownames_to_column("Context") %>%
      pivot_longer(-Context, names_to = "PatientID", values_to = "ZScore") %>%
      left_join(select(p$metadata, PatientID, Group), by = "PatientID") %>%
      mutate(PatientID = factor(PatientID, p$samples),
             Context = factor(Context, rev(rows$Context[order])))
    group_size <- table(p$metadata$Group)
    ggplot(heat, aes(PatientID, Context, fill = ZScore)) +
      geom_raster() +
      facet_grid(. ~ Group, scales = "free_x", space = "free_x",
                 labeller = as_labeller(function(x) paste0(x, " (n = ", group_size[x], ")"))) +
      scale_y_discrete(labels = labels) +
      scale_fill_gradient2(low = "#3266A0", mid = "#FFFFD7", high = "#B62225",
                           name = "Gene z-score") +
      labs(title = paste(p$cohort, "expression of ranked candidate targets"),
           subtitle = paste(gsub("_", " + ", rows$ReferenceRegimen[[1]]),
                            "|", rows$CandidateGridUM[[1]], "µM"),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8),
            axis.text.x = if (length(p$samples) <= 40) element_text(angle = 90, hjust = 1)
                          else element_blank())
  }, height = function() max(360, 120 + 17 * nrow(projection()$rows)))
  output$candidate_cohort_summary <- renderDT({
    datatable(projection()$summary %>%
      select(CellLineName, CandidateRank, CompoundName, TargetGene, CRISPR_Dependency, DependencyClass,
             BlissExpectedIncrementPP, Group, Samples, Detected, MedianExpression, ExpressionUnit),
      rownames = FALSE, options = list(scrollX = TRUE))
  })
  output$candidate_cohort_download <- downloadHandler(
    filename = function() paste0(input$candidate_cohort, "_candidate_expression_summary.csv"),
    content = function(path) write_csv(projection()$summary, path)
  )
}
