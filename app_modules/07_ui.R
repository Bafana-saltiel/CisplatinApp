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
                 DTOutput("ranking")),
        tabPanel("Reference expectation", DTOutput("reference"),
                 h4("Measured dose inputs"), DTOutput("reference_doses")),
        tabPanel("Dose coverage", DTOutput("coverage")),
        tabPanel("Trial registry", p("Registry arms can include other interventions and schedules. Review the linked records for full details."), DTOutput("trials")),
        tabPanel("Target context",
                 p("Annotations for the ten largest expected increments per model, with optional DepMap gene dependency, RNA and mutation data. These annotations describe candidate targets."),
                 DTOutput("target_context")),
        tabPanel("Measured monotherapy",
                 selectInput("mono_drug", "Compound", choices = NULL),
                 plotOutput("measured_plot"), DTOutput("measurements")),
        tabPanel("Fit audit",
                 p("Dose-matched comparison of the original equation and a slope-reoriented diagnostic against Sanger release 8.4 measurements."),
                 DTOutput("audit_summary"), DTOutput("audit_table")),
        tabPanel("Cross-SCC expression", cross_scc_validation_ui())
      )
    )
  )
)
