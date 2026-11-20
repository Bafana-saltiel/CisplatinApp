# Cross-SCC expression sources

The Cross-SCC panel uses SYK, GSK3B, SRC, ERBB2, ERBB3, WEE1, JAK2 and STAT3. The previous revision shipped a different 11-gene panel.

- TCGA-CESC, TCGA-ESCA, TCGA-HNSC and TCGA-LUSC use the UCSC Xena GDC Hub STAR FPKM-UQ matrices at `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-<PROJECT>.star_fpkm-uq.tsv.gz`. Expression is on the source log2(FPKM-UQ + 1) scale. Match the metadata sample barcode exactly when possible, otherwise use a unique barcode with the same participant and sample-type code. TCGA-H7-A6C5 had no corresponding primary tumour in the HNSC matrix and was excluded.
- GSE299125 uses [GEO transcript-level counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299125/suppl/GSE299125_raw_counts.txt.gz). Match versioned transcript IDs to the [GENCODE v44 annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) and sum the mapped transcript-level counts by gene. The resulting gene-level proxy is log2-transformed and standardised within the eligible cohort. This is an exploratory expression measure.
- Sources retrieved 20 November 2026. Cohort metadata and the app's histology filters come from the previous revision.
