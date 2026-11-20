"""Extract cohort RNA values for measured compounds' annotated targets.

Pass the directory containing GENCODE v44, the GSE299125 raw count matrix,
and the four Xena .uq.tsv.gz matrices as the first argument.
"""
import collections
import csv
import gzip
import pathlib
import re
import sys

raw = pathlib.Path(sys.argv[1])
data = pathlib.Path(__file__).resolve().parents[1] / "data"
out = data / "TCGA data"
measured = {row["CompoundID"] for row in csv.DictReader(
    (data / "GDSC2MeasuredViability_Cervix.csv").open())}
targets = {
    gene.strip() for row in csv.DictReader((data / "PortalCompounds.csv").open())
    if row["CompoundID"] in measured
    for gene in re.split("[;,]", row["GeneSymbolOfTargets"])
    if gene.strip()
}

genes, transcripts = {}, {}

# GSE299125 substitutes this reviewed UniProt accession for an Ensembl
# transcript identifier. See https://www.uniprot.org/uniprotkb/Q9H2K2/entry.
geo_aliases = {"Q9H2K2": "TNKS2"}
with gzip.open(raw / "gencode.v44.gtf.gz", "rt") as source:
    for line in source:
        fields = line.split("\t", 8)
        if len(fields) != 9 or fields[2] not in ("gene", "transcript"):
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
        symbol = attrs["gene_name"]
        if symbol not in targets:
            continue
        if fields[2] == "gene":
            genes[attrs["gene_id"].split(".")[0]] = symbol
        else:
            transcripts[attrs["transcript_id"]] = symbol


def write_matrix(cohort, samples, values):
    path = out / f"{cohort}_candidate_expression.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["GeneSymbol", *samples])
        for gene, row in sorted(values.items()):
            writer.writerow([gene, *row])


with gzip.open(raw / "GSE299125_raw_counts.txt.gz", "rt") as source:
    samples = next(source).rstrip().split("\t")[6:]
    totals = [0] * len(samples)
    counts = collections.defaultdict(lambda: [0] * len(samples))
    for line in source:
        row = line.rstrip().split("\t")
        values = list(map(int, row[6:]))
        totals = [a + b for a, b in zip(totals, values)]
        gene = transcripts.get(row[0], geo_aliases.get(row[0]))
        if gene in targets:
            counts[gene] = [a + b for a, b in zip(counts[gene], values)]
    write_matrix("GSE299125", samples, {
        gene: [format(1e6 * value / total, ".10g") for value, total in zip(values, totals)]
        for gene, values in counts.items()
    })

for project in ("CESC", "ESCA", "HNSC", "LUSC"):
    cohort = f"TCGA_{project}"
    with (out / f"{cohort}_sample_metadata.csv").open() as stream:
        metadata = list(csv.DictReader(stream))
    with gzip.open(raw / f"{project}.uq.tsv.gz", "rt") as source:
        header = next(source).rstrip().split("\t")
        index = {sample: i for i, sample in enumerate(header)}
        short = collections.defaultdict(list)
        for sample in header[1:]:
            short[sample[:15]].append(sample)
        samples, indices, values = [], [], {}
        for row in metadata:
            sample = row["PrimaryTumorSampleID"]
            alternatives = short[sample[:15]]
            resolved = sample if sample in index else alternatives[0] if len(alternatives) == 1 else None
            if resolved:
                samples.append(row["PatientID"])
                indices.append(index[resolved])
        for line in source:
            gene = genes.get(line.split("\t", 1)[0].split(".")[0])
            if gene:
                row = line.rstrip("\n").split("\t")
                if gene in values:
                    raise ValueError(f"Duplicate expression rows for {gene} in {cohort}")
                values[gene] = [row[i] for i in indices]
        write_matrix(cohort, samples, values)
