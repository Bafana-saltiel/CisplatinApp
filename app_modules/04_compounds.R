# ============================================================
# 5. AVAILABLE THIRD DRUGS AND ANTI-ANGIOGENIC SUBSET
# ============================================================

Available_Third_Drugs <- bind_rows(
  Cervix_GDSC_All,
  semi_join(GDSC2ResponseCurves, Ovarian_Monotherapy, by = c("ModelID", "CompoundID"))
) %>%
  distinct(CompoundID, CompoundName) %>%
  left_join(
    PortalCompounds %>%
      select(
        CompoundID,
        Synonyms,
        GeneSymbolOfTargets,
        TargetOrMechanism,
        ChEMBLID,
        PubChemCID
      ),
    by = "CompoundID"
  ) %>%
  arrange(CompoundName)

Third_Drug_Choices <- setNames(
  Available_Third_Drugs$CompoundID,
  paste0(
    Available_Third_Drugs$CompoundName,
    " [",
    Available_Third_Drugs$CompoundID,
    "]"
  )
)

# Terms used to identify VEGF/VEGFR and broader angiogenesis-associated
# compounds in the compound catalogue. The final subset is restricted to
# compounds that have an actual GDSC2 response curve in at least one
# supported gynaecological model.
Antiangiogenic_Search_Terms <- c(
  "bevacizumab",
  "avastin",
  "vegf",
  "vegfa",
  "vegfb",
  "vegfc",
  "vegfd",
  "vegfr",
  "vegfr1",
  "vegfr2",
  "vegfr3",
  "flt1",
  "flt4",
  "kdr",
  "vascular endothelial",
  "angiogenesis"
)

make_annotation_text <- function(...) {
  values <- list(...)
  values <- vapply(
    values,
    function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      paste(x, collapse = " ")
    },
    character(1)
  )

  stringr::str_to_lower(paste(values, collapse = " "))
}

identify_matching_terms <- function(search_text, terms) {
  hits <- terms[
    vapply(
      terms,
      function(term) {
        stringr::str_detect(
          search_text,
          stringr::fixed(stringr::str_to_lower(term))
        )
      },
      logical(1)
    )
  ]

  paste(unique(hits), collapse = "; ")
}

Antiangiogenic_Compounds <- Available_Third_Drugs %>%
  rowwise() %>%
  mutate(
    AntiangiogenicSearchText = make_annotation_text(
      CompoundName,
      Synonyms,
      GeneSymbolOfTargets,
      TargetOrMechanism
    ),
    MatchedAntiangiogenicAnnotation = identify_matching_terms(
      AntiangiogenicSearchText,
      Antiangiogenic_Search_Terms
    )
  ) %>%
  ungroup() %>%
  filter(MatchedAntiangiogenicAnnotation != "") %>%
  select(-AntiangiogenicSearchText) %>%
  distinct(CompoundID, .keep_all = TRUE) %>%
  arrange(CompoundName)

Antiangiogenic_Drug_Choices <- setNames(
  Antiangiogenic_Compounds$CompoundID,
  paste0(
    Antiangiogenic_Compounds$CompoundName,
    " [",
    Antiangiogenic_Compounds$CompoundID,
    "]"
  )
)

message(
  "Identified ",
  nrow(Antiangiogenic_Compounds),
  " anti-angiogenic compounds with gynaecological GDSC2 response data."
)

