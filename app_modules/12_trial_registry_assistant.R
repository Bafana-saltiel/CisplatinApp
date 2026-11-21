trial_registry_assistant_ui <- function() {
  tagList(
    h4("Trial registry assistant"),
    actionButton("explain_trial_regimen", "Summarise registry evidence"),
    uiOutput("trial_assistant_answer")
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

trial_registry_assistant_server <- function(input, output, trials_for_regimen) {
  shown_regimen <- eventReactive(input$explain_trial_regimen, input$combo)
  output$trial_assistant_answer <- renderUI({
    req(identical(shown_regimen(), input$combo))
    trials <- summarise_registry_arms(trials_for_regimen(), input$combo)
    validate(need(nrow(trials), "No same-arm registry match was identified for this preset."))
    examples <- trials %>%
      filter(nzchar(AdditionalEntries)) %>%
      arrange(desc(Radiation | AdditionalAgents), NCTId, Arm) %>%
      distinct(NCTId, .keep_all = TRUE) %>%
      slice_head(n = 4L)
    tagList(
      p(sprintf("%s appears in %d curated arm records across %d registered studies.",
                gsub("_", " + ", input$combo, fixed = TRUE), nrow(trials), n_distinct(trials$NCTId))),
      p(sprintf("Radiation is listed in %d arms, additional drugs or biological agents in %d, and other intervention types in %d. These counts can overlap.",
                sum(trials$Radiation), sum(trials$AdditionalAgents), sum(trials$OtherInterventions))),
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
