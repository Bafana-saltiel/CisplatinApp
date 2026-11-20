#!/usr/bin/env python3
"""Retrieve and record example same-arm cisplatin regimens in cervical cancer."""
import csv
import hashlib
import json
import re
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

BASE = "https://clinicaltrials.gov/api/v2"
QUERY = {"query.cond": "Cervical Cancer", "query.intr": "Cisplatin",
         "pageSize": 1000, "countTotal": "true"}
CERVICAL = re.compile(
    r"\b(?:cervical|cervix)\b.*\b(?:cancer|carcinoma|neoplasm|adenocarcinoma|malignancy)\b"
    r"|\b(?:cancer|carcinoma|neoplasm|adenocarcinoma|malignancy)\b.*\b(?:cervical|cervix)\b", re.I)
DRUGS = {
    "Cisplatin": r"cisplatin|platinol|cddp",
    "C225": r"cetuximab|erbitux|c225",
    "Paclitaxel": r"paclitaxel|taxol",
    "Topotecan": r"topotecan|hycamtin",
    "Fluorouracil": r"fluorouracil|5[ -]?fu",
    "Methotrexate": r"methotrexate",
    "Bleomycin": r"bleomycin",
    "Vinorelbine": r"vinorelbine|navelbine",
    "Vinblastine": r"vinblastine",
    "Veliparib": r"veliparib|abt[ -]?888",
}
REGIMENS = [(name,) for name in DRUGS if name not in {"Cisplatin", "Veliparib"}]
REGIMENS.append(("Paclitaxel", "Veliparib"))
# Manually checked registry records excluded from automated matching.
EXCLUDED_STUDIES = {
    "NCT06535607": "Cervical cancer appears only in the monotherapy sub-study.",
    "NCT01405235": "Registered arm label and listed drugs conflict.",
}


def read_json(url):
    with urlopen(url, timeout=100) as response:
        return json.load(response)


def main():
    dirname = Path("data")
    version = read_json(BASE + "/version")
    token = None
    studies = {}
    total = None
    while True:
        query = {**QUERY, **({"pageToken": token} if token else {})}
        payload = read_json(BASE + "/studies?" + urlencode(query))
        if total is None:
            total = payload["totalCount"]
        for study in payload["studies"]:
            record = study["protocolSection"]
            studies[record["identificationModule"]["nctId"]] = record
        token = payload.get("nextPageToken")
        if not token:
            break
    if version != read_json(BASE + "/version") or len(studies) != total:
        raise RuntimeError("Registry changed during download or study retrieval was incomplete.")
    rows = []
    for nct, record in sorted(studies.items()):
        conditions = record.get("conditionsModule", {}).get("conditions", [])
        design = record.get("designModule", {})
        if (nct in EXCLUDED_STUDIES or design.get("studyType") != "INTERVENTIONAL"
                or not any(CERVICAL.search(cond) for cond in conditions)):
            continue
        assigned = record.get("armsInterventionsModule", {})
        arms = assigned.get("armGroups", [])
        if not arms and design.get("designInfo", {}).get("interventionModel") == "SINGLE_GROUP":
            names = [f"{treatment['type']}: {treatment['name']}"
                     for treatment in assigned.get("interventions", [])]
            arms = [{"label": "Single group (registry design)", "interventionNames": names}]
        for arm in arms:
            names = arm.get("interventionNames", [])
            drug_names = " | ".join(name for name in names if name.lower().startswith("drug:"))
            token_is_present = lambda name: re.search(r"(?<![\w-])(?:" + DRUGS[name] +
                                                     r")(?![\w-])", drug_names, re.I)
            if not token_is_present("Cisplatin"):
                continue
            for partners in REGIMENS:
                if all(token_is_present(partner) for partner in partners):
                    rows.append({
                        "Regimen": "Cisplatin_" + "_".join(partners),
                        "NCTId": nct, "Arm": arm["label"],
                        "Interventions": " | ".join(names),
                        "Conditions": " | ".join(conditions),
                        "Title": record["identificationModule"]["briefTitle"],
                        "RecordURL": "https://clinicaltrials.gov/study/" + nct,
                    })
    path = dirname / "clinical-trial-arm-matches.csv"
    if not rows:
        raise RuntimeError("No matching interventional arms were retrieved.")
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {**version, "retrieved": date.today().isoformat(),
                "query_url": BASE + "/studies?" + urlencode(QUERY),
                "filters": "No API study filter. Local selection of INTERVENTIONAL records with cervical disease terms and all regimen drug aliases in the same arm's Drug intervention names. An explicitly single-group study without arm labels uses its study-level interventions.",
                "curated_exclusions": EXCLUDED_STUDIES,
                "returned_studies": len(studies), "matched_arms": len(rows),
                "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (dirname / "clinical-trial-query.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
