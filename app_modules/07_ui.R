ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Cisplatin Clinical Reference Combination Explorer"),
  p("Explore measured GDSC2 monotherapy and calculate Bliss independence expectations for clinically investigated regimens. Experimental combination viability is required to estimate synergy."),
  sidebarLayout(
    sidebarPanel(
      selectInput("combo", "Clinical reference regimen",
                  choices = Clinical_Cisplatin_Combinations_Final$ClinicalCombination,
                  selected = "Cisplatin_Paclitaxel"),
      radioButtons("candidate_mode", "Candidate set",
                   choices = c("Anti-angiogenesis annotations" = "antiangiogenic",
                               "All available compounds" = "all", "Select compounds" = "selected")),
      conditionalPanel("input.candidate_mode === 'selected'",
                       selectizeInput("candidates", "Candidate compounds",
                                      choices = NULL, multiple = TRUE)),
      checkboxGroupInput("doses", "Common candidate dose grid (µM)",
                         choices = c(0.01, 0.03, 0.1, 0.3, 1),
                         selected = c(0.03, 0.1)),
      helpText("Reference drugs use their measured concentrations. Candidates share the selected grid, with log-dose interpolation inside the measured range."),
      checkboxInput("include_hsa", "Include Highest Single Agent (HSA) benchmark", value = FALSE),
      checkboxInput("include_loewe", "Include Loewe additivity expectation", value = FALSE),
      checkboxInput("include_zip", "Include ZIP no-potency-shift expectation", value = FALSE),
      actionButton("run", "Calculate expectations", class = "btn-primary"),
      br(), br(),
      downloadButton("download_ranking", "Download expected increments"),
      br(), br(),
      downloadButton("download_audit", "Download monotherapy fit audit")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Expected increments",
                 p("Within each model on the common grid, the expected increment ranks drugs by their single-agent inhibition after bounding viability to [0, 1]."),
                 uiOutput("hsa_note"),
                 DTOutput("ranking")),
        tabPanel("Top candidates by model",
                 div(style = "overflow-x: auto;",
                     div(style = "min-width: 1000px;",
                         plotOutput("top_candidates_plot", height = "auto")))),
        tabPanel("Reference expectation", DTOutput("reference"),
                 h4("Measured dose inputs"), DTOutput("reference_doses")),
        tabPanel("Dose coverage", DTOutput("coverage")),
        tabPanel("Trial registry", p("Registry arms can include other interventions and schedules. Review the linked records for full details."), DTOutput("trials")),
        tabPanel("Target context",
                 p("Annotations for the ten largest expected increments per model, with optional DepMap gene dependency, RNA and mutation data. These annotations describe candidate targets."),
                 tabsetPanel(
                   tabPanel("CRISPR dependency", plotOutput("target_crispr_plot", height = "auto")),
                   tabPanel("RNA expression", plotOutput("target_expression_plot", height = "auto")),
                   tabPanel("Somatic mutations", plotOutput("target_mutation_plot", height = "auto")),
                   tabPanel("Evidence table", DTOutput("target_context"))
                 )),
        tabPanel("Measured monotherapy",
                 plotOutput("top_monotherapy_plot", height = "auto"),
                 plotOutput("measured_plot", height = "auto"), DTOutput("measurements")),
        tabPanel("Fit audit",
                 p("Dose-matched comparison of the original equation and a slope-reoriented diagnostic against Sanger release 8.4 measurements."),
                 DTOutput("audit_summary"), DTOutput("audit_table")),
        tabPanel("Candidate cohort expression", candidate_cohort_ui()),
        tabPanel("Cross-SCC expression", cross_scc_validation_ui())
      )
    )
  )
)
