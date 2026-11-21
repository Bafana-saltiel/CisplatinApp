# Cisplatin Clinical Reference Combination Explorer

The Shiny app displays measured cervical-line GDSC2 monotherapy and computes Bliss independence expectations for clinically investigated cisplatin regimens. The pharmacological outputs contain measured monotherapy, expected regimen inhibition and expected increments. Experimental combination viability is required to estimate synergy.

## Run and reproduce

From the repository root in R 4.5 or a compatible R installation:

```r
install.packages(c("shiny", "shinythemes", "DT", "dplyr", "stringr", "purrr",
                   "tidyr", "ggplot2", "readr", "tibble", "igraph", "scales"))
shiny::runApp(".")
```

```sh
Rscript scripts/test_measured_screen.R
Rscript scripts/reproduce_screen.R results
python scripts/plot_manuscript.py figures results
```

The batch script reproduces cisplatin–paclitaxel Tables/plots with the 0.03 and 0.10 µM candidate grid and writes the R session information. The plot script requires pandas, NumPy and Matplotlib. The current deployment entry point is `app.R`. Posit Connect Cloud can build from this repository and its `manifest.json`. The [existing deployment](https://bafana-saltiel-cisplatinapp.share.connect.posit.cloud/?utm_source=chatgpt.com) reported Content Unavailable on 20 November 2026. The owner workspace reported exceeded monthly usage credits.

## Measured monotherapy and fit audit

The app ships the cervical subset derived from the [Sanger GDSC2 raw release 8.4](https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/) (dataset dated 24 July 2022, retrieved 20 November 2026). This is a public Sanger release with its own well processing. The reviewer-suggested [DepMap Harmonized GDSC 25Q2 matrix](https://depmap.org/portal/data_page/?tab=allData&releasename=Harmonized%20GDSC%2025Q2&filename=GDSC2Log2ViabilityCollapsedMatrix.csv) was unavailable in the analysis session because the download returned a verification page.

To rebuild the measured table and dose-matched audit with Python 3.11+ and pandas/NumPy:

```sh
curl --fail --location --output GDSC2_public_raw_data_24Jul22.csv \
  'https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC2_public_raw_data_24Jul22.csv'
python scripts/prepare_sanger_viability.py GDSC2_public_raw_data_24Jul22.csv \
  > data/sanger-fit-audit-summary.json
```

The raw input is 2.1 GB. Match Sanger model identifiers and numeric GDSC2 drug identifiers via the archived inputs. Retain single-agent library wells, normalize their signal to the median of at least three NC-1 vehicle wells on the same scan, collapse technical wells and then take the median across plates. The processed file records the plate and well counts and preserves measurements outside [0, 1]. The audit uses pair-matched doses within the archived fit ranges, with 0.1% numerical matching tolerance, and reproduces the original equation and a slope-reoriented diagnostic. It reports equal-weight dose-level residuals. The archived fits supply audit and identifier information. They do not supply viability for the revised screen.

For each backbone drug, the Bliss calculation averages the bounded monotherapy viability over its available tested concentrations. Candidate viability is evaluated on a common absolute concentration grid with log-dose interpolation within measured support. All candidate grid points are required. The mean Bliss independence viability across the Cartesian product is the product of those monotherapy means. Within a model and candidate grid, the expected increment is baseline viability multiplied by one minus candidate viability. This is a deterministic transformation of candidate monotherapy activity. The dose grids represent in vitro scenarios and do not reproduce schedules or doses used in clinical trials.

## Data and provenance

- `data/input-provenance.csv` records input URLs, identifiable releases, SHA-256 digests and earliest matching repository history where available. Historical DepMap snapshots are identified by exact bytes. Their original portal release and acquisition dates are unverified. The reviewer’s named 25Q2 release applies to the suggested measurement matrix.
- The trial panel uses `data/clinical-trial-arm-matches.csv`. `data/clinical-trial-query.json` records the API query, API version, registry data timestamp, retrieval date and local filters. Reproduce a current snapshot with `python scripts/refresh_trial_evidence.py`. Review the arm details for other drugs, alternatives and scheduling. The historical `ctg-studies.csv` is retained for provenance.
- The optional target context tab uses the archived cervical DepMap annotation snapshots. The Cross-SCC tab uses annotated targets from the ranked screen to explore cohort RNA expression.
- The base commit `6b2c0e9128a350285de66081b4da446056c3eaa0` preserves the original implementation for audit. The earlier “observed” combination values came from adding Gaussian noise to Bliss expectations. The revised app removes the synthetic combination scores, the Loewe and ZIP options, and the graph and validation outputs that used synthetic interaction labels.

### Candidate cohort expression

After calculating expectations, open **Candidate cohort expression**. Choose the cohort and the number of ranked candidates per cell line. The heatmap and export use the target rows from the current run’s Target context tab. Rank 1 is the default and projects the top-ranked compound from each cell line, including in Select compounds mode. The optional single-target view selects the highest available CRISPR dependency within each compound–model pair and uses alphabetical order for ties or when no dependency score exists. GSE299125 can be grouped by HIV status among invasive tumours or by HIV and disease status across the sampled cervix groups. TCGA primary tumours can be grouped by broad or individual stage, SCC status, histologic diagnosis, grade, age at diagnosis, recorded sex, anatomical site or recorded vital status. The menu shows options with at least two recorded groups in the selected cohort. Missing values are excluded from the selected grouping. The heatmap clusters samples within their selected groups and the summary can be exported as CSV. See [expression sources](data/cross-scc-expression-sources.md).

### Cross-SCC Reactome pathways

Calculate expectations, choose a Reactome pathway in **Cross-SCC expression**, and run the RNA analysis. The search lists pathways overlapping the current ranked targets. The same measured target intersection is used across all five cohorts and every result view. **All screened targets** provides a broader view. Programme scores average the within-cohort RNA z-scores of selected targets mapped to cell cycle, PI3K–AKT, cell death, DNA repair and replication, immunity, metabolism, stress, RTK and VEGF signalling, and the extracellular matrix. The bundled [Reactome v98 gene sets](https://reactome.org/download-data) use stable human pathway identifiers. GSE299125 candidate CPM values are log-transformed before standardisation.

### Ovarian models and experimental combination uploads

Choose **Ovarian/fallopian tube** under **Cancer type** to screen that DepMap lineage using measured Sanger GDSC2 monotherapy. The clinical registry evidence and fit audit are from the cervical cancer analysis. Ovarian target rows use available annotations, with missing model-specific DepMap evidence shown as unavailable.

Use **Experimental combinations** to view the NCI ALMANAC ovarian-panel example or upload your own two-drug dose matrix. The example is downloadable in the same format the upload accepts. Each CSV row supplies a model, experiment, two drug names, both concentrations in µM, matched single-agent viability fractions at those doses and the measured combination viability fraction. Use a vehicle-relative scale (1 is vehicle control). Each positive-dose combination row is compared to the Bliss and HSA expectations in percentage points. This tab analyses the experimental pair independently of the monotherapy-based third-drug screen.

The example uses ALMANAC test signal divided by matched control signal. The single-agent values come from the same study, screener, test date, plate, drug sample and dose. To rebuild the example, pass the [NCI ALMANAC ZIP](https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip) and `data/ALMANAC_ovarian_cisplatin_paclitaxel.csv` to `scripts/prepare_almanac_benchmark.py`. To rebuild ovarian monotherapy, run `scripts/prepare_sanger_viability.py` with the Sanger raw CSV and `--lineage Ovary`.

The ovarian benchmark also compares expected inhibition computed from GDSC2 monotherapy with observed NCI ALMANAC inhibition for matching models and dose pairs. GDSC2 responses are interpolated on log concentration within the measured range. The CSV includes dose pairs excluded for missing model matches or dose coverage. Comparisons use independent assays and may reflect assay differences.
