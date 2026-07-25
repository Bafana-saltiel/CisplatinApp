# ============================================================
# 2. CISPLATIN REFERENCE COMBINATIONS
# ============================================================

Clinical_Cisplatin_Combinations_Final <- tibble::tribble(
  ~ClinicalCombination, ~TherapeuticClass, ~DrugNames, ~CompoundIDs,
  "Cisplatin_C225", "EGFR-targeted antibody", c("CISPLATIN", "C225"), c("DPC-001793", "DPC-001418"),
  "Cisplatin_Paclitaxel", "Microtubule inhibitor", c("CISPLATIN", "PACLITAXEL"), c("DPC-001793", "DPC-004880"),
  "Cisplatin_Topotecan", "Topoisomerase I inhibitor", c("CISPLATIN", "TOPOTECAN"), c("DPC-001793", "DPC-006564"),
  "Cisplatin_Fluorouracil", "Antimetabolite", c("CISPLATIN", "FLUOROURACIL"), c("DPC-001793", "DPC-002828"),
  "Cisplatin_Methotrexate", "Antifolate", c("CISPLATIN", "METHOTREXATE"), c("DPC-001793", "DPC-000234"),
  "Cisplatin_Bleomycin", "Cytotoxic chemotherapy", c("CISPLATIN", "BLEOMYCIN"), c("DPC-001793", "DPC-001132"),
  "Cisplatin_Vinorelbine", "Microtubule inhibitor", c("CISPLATIN", "VINORELBINE"), c("DPC-001793", "DPC-002684"),
  "Cisplatin_Vinblastine", "Vinca alkaloid", c("CISPLATIN", "VINBLASTINE"), c("DPC-001793", "DPC-006819"),
  "Cisplatin_Paclitaxel_Veliparib", "PARP inhibitor combination", c("CISPLATIN", "PACLITAXEL", "VELIPARIB"), c("DPC-001793", "DPC-004880", "DPC-000236")
)

# ============================================================
# 3. PREPARE CERVIX DATA
# ============================================================

Cervix_Model <- Model %>%
  filter(OncotreeLineage == "Cervix") %>%
  select(ModelID, CellLineName, StrippedCellLineName, DepmapModelType)

Cervix_GDSC_All <- GDSC2ResponseCurves %>%
  inner_join(Cervix_Model, by = "ModelID") %>%
  distinct(
    CompoundID, CompoundName, SampleID, ModelID, CellLineName,
    EC50, LowerAsymptote, UpperAsymptote, Slope,
    MinimumDose, MaximumDose, DoseUnit,
    .keep_all = TRUE
  )

