# ============================================================
# 7. SHINY UI
# ============================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Cisplatin Clinical Reference Combination Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput(
        "combo",
        "Select clinical cisplatin reference combination:",
        choices = Clinical_Cisplatin_Combinations_Final$ClinicalCombination,
        selected = "Cisplatin_Paclitaxel"
      ),

      hr(),

      radioButtons(
        "third_drug_mode",
        "Third-drug screening mode:",
        choices = c(
          "Screen all cervical GDSC2 compounds" = "all",
          "Screen all anti-angiogenic cervical compounds" = "antiangiogenic",
          "Select third drugs manually" = "selected"
        ),
        selected = "antiangiogenic"
      ),

      conditionalPanel(
        condition = "input.third_drug_mode == 'antiangiogenic'",
        wellPanel(
          tags$b("Anti-angiogenic screening set"),
          tags$p(
            paste0(
              nrow(Cervical_Antiangiogenic_Compounds),
              " compounds with VEGF-, VEGFR-, KDR-, FLT1-, FLT4- or angiogenesis-related annotations and cervical GDSC2 response data will be selected automatically."
            )
          ),
          tags$p(
            "The subset is generated dynamically from PortalCompounds.csv and Cervix_GDSC_All, so only compounds tested in at least one cervical cancer model are included."
          )
        )
      ),

      conditionalPanel(
        condition = "input.third_drug_mode == 'selected'",
        selectizeInput(
          "selected_third_drugs",
          "Search and select third drugs:",
          choices = Third_Drug_Choices,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Type drug name or CompoundID...",
            maxOptions = 1000,
            plugins = list("remove_button")
          )
        ),
        helpText(
          "Each selected drug is evaluated independently as monotherapy and as an addition to the selected clinical reference backbone."
        )
      ),

      hr(),

      selectInput(
        "synergy_method",
        "Synergy model:",
        choices = c(
          "Bliss independence" = "bliss",
          "Highest Single Agent (HSA)" = "hsa",
          "Loewe additivity" = "loewe",
          "ZIP approximation" = "zip"
        ),
        selected = "bliss"
      ),

      helpText(
        paste(
          "The selected model defines the expected non-interacting response.",
          "Synergy score = observed response - model-specific expected response.",
          "The ZIP option is a transparent curve-based approximation because",
          "the app reconstructs response surfaces from fitted single-agent curves."
        )
      ),

      hr(),

      checkboxInput(
        "enable_gnn",
        "Run graph neural network drug prioritisation",
        value = FALSE
      ),

      conditionalPanel(
        condition = "input.enable_gnn == true",
        selectInput(
          "graph_engine",
          "Graph-learning engine:",
          choices = c(
            "Automatic: GCN when torch is available" = "auto",
            "Two-layer graph convolutional network (requires torch)" = "gcn",
            "Transparent graph-propagation baseline" = "propagation"
          ),
          selected = "auto"
        ),
        sliderInput(
          "graph_edge_threshold",
          "Minimum graph edge weight:",
          min = 0.05, max = 0.80, value = 0.20, step = 0.05
        ),
        numericInput(
          "gnn_epochs",
          "GCN training epochs:",
          value = 250, min = 50, max = 2000, step = 50
        ),
        numericInput(
          "gnn_hidden_dim",
          "GCN hidden units:",
          value = 16, min = 4, max = 128, step = 4
        ),
        helpText(
          paste0(
            "Each candidate drug is represented as a graph node. Edges combine shared targets, ",
            "shared mechanism/pathway terms, CRISPR similarity, expression similarity and mutation context. ",
            "The GCN predicts whether a candidate exceeds the model-specific clinical reference. ",
            "R torch detected: ", has_torch, "."
          )
        )
      ),

      hr(),

      checkboxInput(
        "enable_validation_agent",
        "Run Validation Agent for graph model",
        value = FALSE
      ),

      conditionalPanel(
        condition = "input.enable_validation_agent == true",
        numericInput(
          "validation_folds",
          "Stratified cross-validation folds:",
          value = 5, min = 2, max = 10, step = 1
        ),
        sliderInput(
          "validation_threshold",
          "Classification decision threshold:",
          min = 0.10, max = 0.90, value = 0.50, step = 0.05
        ),
        helpText(
          paste(
            "The Validation Agent performs stratified node-level cross-validation",
            "within each cervical model and benchmarks the graph model against",
            "a pharmacology-only logistic model and a prevalence baseline."
          )
        )
      ),

      hr(),

      checkboxGroupInput(
        "doses",
        "Dose levels:",
        choices = c(0.001, 0.01, 0.1, 1),
        selected = c(0.001, 0.01, 0.1, 1)
      ),

      numericInput("noise", "Noise SD:", value = 5, min = 0, max = 20, step = 1),
      actionButton("run", "Run analysis", class = "btn-primary"),

      br(), br(),
      downloadButton(
        "download_antiangiogenic_set",
        "Download anti-angiogenic compound set"
      ),
      br(), br(),
      downloadButton("download_reference", "Download reference summary"),
      br(), br(),
      downloadButton("download_monotherapy", "Download monotherapy summary"),
      br(), br(),
      downloadButton("download_monotherapy_dose", "Download monotherapy dose-response"),
      br(), br(),
      downloadButton("download_ranking", "Download third-drug ranking"),
      br(), br(),
      downloadButton("download_top10", "Download top 10"),
      br(), br(),
      downloadButton("download_functional", "Download functional CRISPR"),
      br(), br(),
      downloadButton("download_molecular", "Download molecular expression"),
      br(), br(),
      downloadButton("download_biomarker", "Download optional biomarker"),
      br(), br(),
      conditionalPanel(
        condition = "input.enable_gnn == true",
        downloadButton("download_graph_ranking", "Download graph ranking"),
        br(), br(),
        downloadButton("download_graph_edges", "Download graph edges"),
        br(), br(),
        downloadButton("download_graph_diagnostics", "Download graph diagnostics"),
        br(), br(),
        conditionalPanel(
          condition = "input.enable_validation_agent == true",
          downloadButton("download_validation_summary", "Download validation summary"),
          br(), br(),
          downloadButton("download_validation_predictions", "Download CV predictions"),
          br(), br(),
          downloadButton("download_validation_folds", "Download fold metrics"),
          br(), br(),
          downloadButton("download_decision_curve", "Download decision-curve data")
        )
      )
    ),

    mainPanel(
      verbatimTextOutput("run_status"),

      tabsetPanel(
        tabPanel("Eligible models", DTOutput("eligible_table")),
        tabPanel("Reference summary", DTOutput("reference_table")),
        tabPanel("Selected third drugs", DTOutput("third_drug_selection_table")),

        tabPanel(
          "Monotherapy response",
          h4("How strongly does each third drug inhibit the matched cervical cancer model by itself?"),
          DTOutput("monotherapy_table"),
          plotOutput("monotherapy_plot", height = "850px")
        ),

        tabPanel(
          "Monotherapy dose-response",
          h4("Predicted single-drug inhibition across the selected dose levels"),
          DTOutput("monotherapy_dose_table"),
          plotOutput("monotherapy_dose_plot", height = "850px")
        ),

        tabPanel("Third-drug ranking", DTOutput("ranking_table")),
        tabPanel("Top candidates per model", DTOutput("top10_table")),
        tabPanel("Synergy plot", plotOutput("top10_plot", height = "850px")),

        tabPanel(
          "Mono vs combination",
          h4("Comparison of intrinsic monotherapy response with observed combination response"),
          plotOutput("mono_combination_plot", height = "850px")
        ),

        tabPanel(
          "Graph ML ranking",
          h4("Graph-based prediction of the best third drug"),
          p("Candidate drugs are nodes; edges encode shared targets, mechanisms, dependencies, expression and mutation context."),
          DTOutput("graph_ranking_table"),
          plotOutput("graph_ranking_plot", height = "850px")
        ),

        tabPanel(
          "Drug graph",
          selectInput("graph_model_view", "Cervical model to visualise:", choices = NULL),
          plotOutput("drug_graph_plot", height = "850px"),
          DTOutput("graph_diagnostics_table")
        ),

        tabPanel(
          "Validation Agent",
          h4("Graph-model validation and benchmarking"),
          p("Reports accuracy, sensitivity, specificity, ROC-AUC, F1-score, MCC, confusion matrices, decision-curve analysis, clustering validity, and stratified cross-validation."),
          DTOutput("validation_summary_table"),
          h4("Cross-validation metrics by cervical model and fold"),
          DTOutput("validation_fold_table")
        ),

        tabPanel(
          "ROC and confusion matrix",
          plotOutput("validation_roc_plot", height = "650px"),
          plotOutput("validation_confusion_plot", height = "650px")
        ),

        tabPanel(
          "Decision Curve Analysis",
          plotOutput("decision_curve_plot", height = "700px")
        ),

        tabPanel(
          "Graph clustering validation",
          h4("Agreement between graph-derived clusters and pharmacological interaction labels"),
          DTOutput("clustering_validation_table"),
          plotOutput("clustering_validation_plot", height = "650px")
        ),

        tabPanel(
          "Functional: CRISPR dependency",
          h4("Does the model require the drug target for survival?"),
          DTOutput("functional_table"),
          plotOutput("functional_plot", height = "750px")
        ),

        tabPanel(
          "Molecular: target expression",
          h4("Is the target actually expressed in that model?"),
          DTOutput("molecular_table"),
          plotOutput("molecular_plot", height = "750px")
        ),

        tabPanel(
          "Optional biomarker: mutations",
          h4("Is dependency or synergy improvement associated with a genomic alteration?"),
          DTOutput("biomarker_table"),
          plotOutput("biomarker_plot", height = "650px")
        )
      )
    )
  )
)
