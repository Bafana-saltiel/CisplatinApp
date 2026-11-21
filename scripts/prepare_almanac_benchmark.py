"""Extract same-plate ovarian cisplatin–paclitaxel measurements from the NCI ALMANAC ZIP."""
import argparse
import csv
import io
import statistics
import zipfile
from collections import defaultdict

CISPLATIN, PACLITAXEL = "119875", "125973"
DRUGS = {CISPLATIN, PACLITAXEL}


def extract(archive):
    controls = defaultdict(list)
    pairs = []
    with zipfile.ZipFile(archive) as zipped:
        with zipped.open("ComboDrugGrowth_Nov2017.csv") as raw:
            for row in csv.DictReader(io.TextIOWrapper(raw)):
                first, second = row["NSC1"], row["NSC2"]
                if row["PANEL"] != "Ovarian Cancer" or row["VALID"] != "Y" or first not in DRUGS:
                    continue
                if second and {first, second} != DRUGS:
                    continue
                if not row["TESTVALUE"] or not row["CONTROLVALUE"] or float(row["CONTROLVALUE"]) <= 0:
                    continue
                key = tuple(row[name] for name in ("CELLNAME", "SCREENER", "STUDY", "TESTDATE", "PLATE"))
                viability = float(row["TESTVALUE"]) / float(row["CONTROLVALUE"])
                if row["CONCUNIT1"] != "M" or (second and row["CONCUNIT2"] != "M"):
                    raise ValueError("Expected molar concentrations")
                dose = float(row["CONC1"])
                if not second:
                    controls[key, first, row["SAMPLE1"], dose].append(viability)
                    continue
                pair = {first: (row["SAMPLE1"], dose), second: (row["SAMPLE2"], float(row["CONC2"]))}
                pairs.append((key, pair, viability, row["COMBODRUGSEQ"]))

    if not pairs:
        raise ValueError("No valid ovarian cisplatin–paclitaxel measurements found")
    yield ("Model", "Experiment", "DrugA", "DrugB", "DoseA_uM", "DoseB_uM",
           "ViabilityA", "ViabilityB", "CombinationViability", "SourceRow")
    for key, pair, viability, source_row in pairs:
        sample_a, dose_a = pair[CISPLATIN]
        sample_b, dose_b = pair[PACLITAXEL]
        a = controls[key, CISPLATIN, sample_a, dose_a]
        b = controls[key, PACLITAXEL, sample_b, dose_b]
        if not a or not b:
            raise ValueError(f"Missing matched single-agent controls for record {source_row}")
        yield (key[0], "|".join(key[1:]), "Cisplatin", "Paclitaxel", dose_a * 1e6, dose_b * 1e6,
               statistics.median(a), statistics.median(b), viability, source_row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    parser.add_argument("output")
    args = parser.parse_args()
    rows = list(extract(args.archive))
    with open(args.output, "w", newline="") as output:
        csv.writer(output, lineterminator="\n").writerows(rows)
