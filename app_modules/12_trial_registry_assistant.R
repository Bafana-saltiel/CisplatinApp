trial_registry_assistant_ui <- function() {
  tagList(
    h4("Trial registry assistant"),
    actionButton("explain_trial_regimen", "Summarise registry evidence"),
    uiOutput("trial_assistant_answer"),
    tags$button(id = "explain_trial_qwen", type = "button", class = "btn btn-default",
                disabled = NA, "Explain with Qwen"),
    tags$div(id = "trial_qwen_facts", hidden = NA,
             `data-facts` = jsonlite::toJSON(trial_registry_facts, auto_unbox = TRUE)),
    p("Optional. Qwen downloads on first use and runs in your WebGPU browser. The registry summary above works without it."),
    tags$div(
      id = "trial_qwen_progress_panel", hidden = NA,
      tags$label("Qwen model download and loading", `for` = "trial_qwen_progress"),
      tags$progress(id = "trial_qwen_progress", max = 100,
                    style = "display: block; width: 100%; max-width: 28rem; height: 1.25rem;"),
      tags$p(id = "trial_qwen_progress_label")
    ),
    tags$p("Preparing the Qwen controls…", id = "trial_qwen_status", role = "status"),
    tags$div(id = "trial_qwen_answer", style = "white-space: pre-wrap;"),
    tags$script(src = "trial-registry-qwen.js?v=4", defer = NA)
  )
}

summarise_registry_arms <- function(trials, regimen) {
  aliases <- c(Cisplatin = "cisplatin|platinol|cddp",
               Paclitaxel = "paclitaxel|taxol",
               Topotecan = "topotecan|hycamtin",
               Fluorouracil = "fluorouracil|5[ -]?fu",
               Vinorelbine = "vinorelbine|navelbine",
               Veliparib = "veliparib|abt[ -]?888")
  components <- strsplit(regimen, "_", fixed = TRUE)[[1]]
  patterns <- unname(aliases[components])
  patterns[is.na(patterns)] <- components[is.na(patterns)]
  regimen_pattern <- stringr::regex(paste0("\\b(?:", paste(patterns, collapse = "|"), ")\\b"),
                                   ignore_case = TRUE)
  trials %>%
    distinct(NCTId, Arm, Interventions) %>%
    mutate(ArmRow = row_number()) %>%
    tidyr::separate_longer_delim(Interventions, delim = " | ") %>%
    mutate(
      Category = toupper(trimws(sub(":.*$", "", Interventions))),
      Intervention = trimws(sub("^[^:]*:", "", Interventions)),
      ReferenceDrug = Category == "DRUG" & stringr::str_detect(Intervention, regimen_pattern)
    ) %>%
    group_by(ArmRow, NCTId, Arm) %>%
    summarise(
      Radiation = any(Category == "RADIATION", na.rm = TRUE),
      AdditionalAgents = any(Category %in% c("DRUG", "BIOLOGICAL") & !ReferenceDrug, na.rm = TRUE),
      OtherInterventions = any(!(Category %in% c("DRUG", "BIOLOGICAL", "RADIATION")) &
                                 nzchar(Intervention), na.rm = TRUE),
      AdditionalEntries = paste(Interventions[which(!ReferenceDrug)], collapse = " | "),
      .groups = "drop"
    )
}

registry_evidence <- function(trials, regimen) {
  arms <- summarise_registry_arms(trials, regimen)
  examples <- arms %>%
    filter(nzchar(AdditionalEntries)) %>%
    arrange(desc(Radiation | AdditionalAgents), NCTId, Arm) %>%
    distinct(NCTId, .keep_all = TRUE) %>%
    slice_head(n = 4L) %>%
    select(NCTId, Arm, AdditionalEntries)
  list(
    regimen = gsub("_", " + ", regimen, fixed = TRUE),
    arms = nrow(arms), studies = n_distinct(arms$NCTId),
    radiation_arms = sum(arms$Radiation), additional_agent_arms = sum(arms$AdditionalAgents),
    other_intervention_arms = sum(arms$OtherInterventions),
    examples = examples
  )
}

registry_regimens <- Clinical_Cisplatin_Combinations_Final$ClinicalCombination
trial_registry_facts <- setNames(lapply(registry_regimens, function(regimen) {
  registry_evidence(filter(Trial_Arm_Matches, Regimen == regimen), regimen)
}), registry_regimens)

trial_registry_assistant_server <- function(input, output) {
  shown_regimen <- eventReactive(input$explain_trial_regimen, input$combo)
  output$trial_assistant_answer <- renderUI({
    req(identical(shown_regimen(), input$combo))
    facts <- trial_registry_facts[[input$combo]]
    validate(need(facts$arms > 0, "No same-arm registry match was identified for this preset."))
    examples <- facts$examples
    tagList(
      p(sprintf("%s appears in %d curated arm records across %d registered studies.",
                facts$regimen, facts$arms, facts$studies)),
      p(sprintf("Radiation is listed in %d arms, additional drugs or biological agents in %d, and other intervention types in %d. These counts can overlap.",
                facts$radiation_arms, facts$additional_agent_arms, facts$other_intervention_arms)),
      if (nrow(examples)) tagList(
        h5("Examples of additional arm entries"),
        tags$ul(lapply(seq_len(nrow(examples)), function(i) {
          id <- examples$NCTId[[i]]
          tags$li(tags$a(id, href = paste0("https://clinicaltrials.gov/study/", id)),
                  paste(examples$Arm[[i]], examples$AdditionalEntries[[i]], sep = " — "))
        }))
      ),
      p("Registry arm matches describe what was studied. Assessments and procedures are shown as registered, and the summary does not assess efficacy.")
    )
  })
}

qwen_panel_ui <- function(id) {
  tags$details(
    tags$summary("Qwen explanation"),
    uiOutput(paste0("qwen_context_", id))
  )
}

qwen_panel_server <- function(output, id, title, get_rows, guidance, aggregate_only = FALSE) {
  output[[paste0("qwen_context_", id)]] <- renderUI({
    rows <- get_rows()
    validate(need(nrow(rows), "There are no results to explain for this selection."))
    numeric_rows <- select(rows, where(is.numeric))
    ranges <- lapply(numeric_rows, function(values) {
      values <- values[is.finite(values)]
      if (length(values)) round(c(min = min(values), mean = mean(values), max = max(values)), 3) else NULL
    })
    ranges <- Filter(Negate(is.null), ranges)
    facts <- list(
      panel = title, guidance = guidance, row_count = nrow(rows),
      numeric_ranges = lapply(ranges, as.list),
      example_rows = if (!aggregate_only) head(select(rows, where(~ is.atomic(.x))), 5L) else NULL
    )
    tags$div(
      class = "qwen-panel", `data-context` = jsonlite::toJSON(facts, auto_unbox = TRUE, na = "null"),
      tags$button(type = "button", class = "btn btn-default qwen-panel-button", "Explain with Qwen"),
      tags$p(class = "qwen-panel-status", role = "status"),
      tags$p(class = "qwen-panel-answer", style = "white-space: pre-wrap;")
    )
  })
}
