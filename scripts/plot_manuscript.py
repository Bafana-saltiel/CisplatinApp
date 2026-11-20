#!/usr/bin/env python3
"""Create figures for the measured-monotherapy manuscript."""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

output = Path(sys.argv[1])
analysis = Path(sys.argv[2])
output.mkdir(exist_ok=True, parents=True)
audit = pd.read_csv("data/GDSC2MonotherapyFitAudit_Cervix.csv")
stats = json.loads(Path("data/sanger-fit-audit-summary.json").read_text())
samples = pd.read_csv(analysis / "seven_candidates.csv")

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), layout="constrained", sharey=True)
limits = (min(-0.05, audit.ViabilityFraction.min()), max(1.05, audit.ViabilityFraction.max()))
bins = []
for ax, key, name in zip(axes, ["OriginalFit", "SlopeReorientedFit"],
                         ["A   Original equation", "B   Reoriented slope diagnostic"]):
    dots = ax.hexbin(audit.ViabilityFraction, audit[key], gridsize=48,
                     mincnt=1, bins="log", cmap="viridis")
    bins.append(dots)
    ax.axline((0, 0), slope=1, color="#933431", linestyle="--", linewidth=1)
    ax.set(xlim=limits, ylim=(-.05, 1.08), title=name, xlabel="Measured viability fraction")
    ax.text(.98, .96, f"MAE {stats[key]['mae']:.3f}", ha="right", va="top",
            transform=ax.transAxes, bbox={"facecolor": "white", "alpha": .9, "edgecolor": "none"})
axes[0].set_ylabel("Reconstructed viability fraction")
max_count = max(artist.get_array().max() for artist in bins)
for artist in bins:
    artist.set_clim(1, max_count)
fig.colorbar(bins[0], ax=axes, label="Dose records")
fig.savefig(output / "monotherapy_fit_audit.png", dpi=300)
plt.close(fig)

names = ["DACTINOMYCIN", "ABT737", "XAV-939", "FULVESTRANT", "BMS-754807"]
fig, ax = plt.subplots(figsize=(8, 2.9), layout="constrained")
for index, name in enumerate(names):
    values = np.sort(samples.loc[samples.CompoundName.eq(name), "BlissExpectedIncrementPP"])
    jitter = np.linspace(-.13, .13, len(values))
    ax.scatter(values, index + jitter, s=23, color="#007b81", alpha=.72, linewidths=0)
    ax.scatter(np.median(values), index, s=90, color="#133451", marker="|", linewidths=2)
ax.set(yticks=range(len(names)), yticklabels=[n.title() if n in {"DACTINOMYCIN", "FULVESTRANT"} else n for n in names],
       xlabel="Expected inhibition increment (percentage points)", xlim=(-2, 100))
ax.invert_yaxis()
ax.grid(axis="x", alpha=.2)
fig.savefig(output / "candidate_expected_increments.png", dpi=300)
