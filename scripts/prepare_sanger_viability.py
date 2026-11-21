#!/usr/bin/env python3
"""Derive vehicle-relative monotherapy viability from Sanger GDSC2 release 8.4."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SANGER_RELEASE = "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/"
RAW_FILE = "GDSC2_public_raw_data_24Jul22.csv"
RAW_BYTES = 2112924313


def prepare(path, data_dir, lineage="Cervix"):
    models = pd.read_csv(data_dir / "Model.csv")
    models = models.loc[models.OncotreeLineage.eq({"Cervix": "Cervix", "Ovary": "Ovary/Fallopian Tube"}[lineage]), ["SangerModelID", "ModelID"]].dropna()
    curves = pd.read_csv(data_dir / "GDSC2ResponseCurves.csv")
    curves = curves[curves.ModelID.isin(models.ModelID)].copy()
    curves["DRUG_ID"] = curves.SampleID.str.extract(r"^GDSC2:(\d+)$").astype(int)
    assert not curves.duplicated(["ModelID", "DRUG_ID"]).any()
    cols = ["SANGER_MODEL_ID", "BARCODE", "SCAN_ID", "DURATION", "ASSAY", "POSITION",
            "TAG", "DRUG_ID", "CONC", "INTENSITY"]
    chunks = pd.read_csv(path, usecols=cols, chunksize=500_000)
    wells = pd.concat(chunk[chunk.SANGER_MODEL_ID.isin(models.SangerModelID)] for chunk in chunks)
    wells = wells.dropna(subset=["INTENSITY", "BARCODE", "SCAN_ID", "POSITION"])
    plate = ["SANGER_MODEL_ID", "BARCODE", "SCAN_ID"]
    well = plate + ["POSITION"]
    vehicle = wells[wells.TAG.eq("NC-1")].drop_duplicates(well)
    controls = vehicle.groupby(plate).INTENSITY.agg(VehicleMedian="median", VehicleWells="size")
    controls = controls[(controls.VehicleMedian > 0) & (controls.VehicleWells >= 3)]
    # A library single-agent well cannot have another drug treatment tag.
    treated = wells[wells.TAG.str.fullmatch(r"L\d+-D\d+-S", na=False)].copy()
    drug_tags = wells[wells.DRUG_ID.notna() & wells.CONC.gt(0)].groupby(well).DRUG_ID.nunique()
    treated = treated.merge(drug_tags[drug_tags.eq(1)].rename("SingleDrug"), on=well)
    treated = treated[treated.CONC.gt(0)].drop_duplicates(well)
    treated = treated.merge(controls, on=plate, validate="many_to_one")
    treated = treated.merge(models, left_on="SANGER_MODEL_ID", right_on="SangerModelID",
                            validate="many_to_one")
    treated = treated.merge(curves, on=["ModelID", "DRUG_ID"], validate="many_to_one")
    treated["ViabilityFraction"] = treated.INTENSITY / treated.VehicleMedian
    treated = treated[np.isfinite(treated.ViabilityFraction)]
    treated["CONC"] = treated.CONC.map(lambda dose: float(f"{dose:.10g}"))
    keys = ["ModelID", "CompoundID", "SampleID", "CONC"]
    # Equal weight per plate after collapsing technical wells on each plate.
    per_plate = treated.groupby(keys + ["BARCODE", "SCAN_ID"]).ViabilityFraction.agg(
        PlateViability="median", NWells="size").reset_index()
    observed = per_plate.groupby(keys).agg(
        ViabilityFraction=("PlateViability", "median"),
        NPlates=("PlateViability", "size"),
        NWells=("NWells", "sum"),
    ).reset_index().rename(columns={"CONC": "DoseUM"})
    if observed.empty:
        raise ValueError("No Sanger measurements matched the selected lineage and GDSC2 drug IDs")
    observed.to_csv(data_dir / f"GDSC2MeasuredViability_{lineage}.csv", index=False, float_format="%.12g")
    # The original equation is reproduced verbatim for the reviewer-requested audit.
    audit = observed.merge(curves, on=["ModelID", "CompoundID", "SampleID"], validate="many_to_one")
    audit = audit[(audit.DoseUnit == "uM") &
                  audit.DoseUM.between(audit.MinimumDose * 0.999, audit.MaximumDose * 1.001)].copy()
    with np.errstate(over="ignore"):
        log_ratio = np.log(audit.DoseUM) - np.log(audit.EC50)
        for label, sign in [("OriginalFit", 1), ("SlopeReorientedFit", -1)]:
            predicted = audit.LowerAsymptote + (audit.UpperAsymptote - audit.LowerAsymptote) / (
                1 + np.exp(sign * audit.Slope * log_ratio))
            audit[label] = predicted
    audit_columns = ["ModelID", "CompoundID", "SampleID", "CompoundName", "DoseUM",
                     "ViabilityFraction", "NPlates", "NWells", "OriginalFit", "SlopeReorientedFit"]
    audit[audit_columns].to_csv(data_dir / f"GDSC2MonotherapyFitAudit_{lineage}.csv", index=False,
                               float_format="%.12g")
    summary = {"matched_doses": len(audit), "models": int(audit.ModelID.nunique()),
               "model_compound_pairs": len(audit[["ModelID", "CompoundID"]].drop_duplicates()),
               "measurement_rows": len(observed), "matched_wells": len(treated),
               "durations_days": sorted(treated.DURATION.dropna().unique().tolist()),
               "assays": sorted(treated.ASSAY.dropna().unique().tolist())}
    for label in ("OriginalFit", "SlopeReorientedFit"):
        error = audit[label] - audit.ViabilityFraction
        summary[label] = {"mae": float(error.abs().mean()), "rmse": float(np.sqrt((error**2).mean())),
                          "bias": float(error.mean())}
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_file", type=Path)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--lineage", choices=("Cervix", "Ovary"), default="Cervix")
    args = parser.parse_args()
    if args.raw_file.stat().st_size != RAW_BYTES:
        parser.error(f"Expected the complete {RAW_FILE} ({RAW_BYTES} bytes)")
    with args.raw_file.open("rb") as source:
        sha256 = hashlib.file_digest(source, "sha256").hexdigest()
    summary = prepare(args.raw_file, args.data_dir, args.lineage)
    summary.update(source=SANGER_RELEASE + RAW_FILE, raw_sha256=sha256)
    print(json.dumps(summary, indent=2, allow_nan=False))
