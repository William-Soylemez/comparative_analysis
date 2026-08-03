#!/usr/bin/env python3
"""
Overlay comparison of three IsoRank alignment score distributions on one
chart: the paper's own bakers<->rat reference alignment (paper_reference_check/,
using whichever similarity preprocessing the paper's own rat-bakers.tsv
already encodes) against our ncrassa<->calbicans alignment run twice, once
per similarity metric (E-score = -log10(evalue), bitscore = raw bitscore).

log10(score) on x (all three use the same convention: higher/closer to 0 =
stronger match), density (not raw count) on y so the three curves are
comparable despite different n (rat: 6478 pairs, ours: 4555 pairs each) --
each curve integrates to 1 over its own histogram.

Usage:
    python3 plot_rat_vs_ours_overlay.py
Expects paper_reference_check/bakers_rat_alignment_scored.tsv,
        ncrassa_calbicans_alignment_scored.tsv,
        ncrassa_bs_calbicans_bs_alignment_scored.tsv
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2a78d6", "#1baf7a", "#eda100"]
BIN_WIDTH = 0.5

series = [
    ("rat (paper reference, bakers<->rat)", "paper_reference_check/bakers_rat_alignment_scored.tsv"),
    ("ncrassa<->calbicans, E-score (-log10 evalue)", "results/ncrassa_calbicans_alignment_scored.tsv"),
    ("ncrassa<->calbicans, bitscore", "results/ncrassa_bs_calbicans_bs_alignment_scored.tsv"),
]

data = []
for label, path in series:
    scores = pd.read_csv(path, sep="\t")["score"].to_numpy()
    data.append((label, np.log10(scores)))

lo = min(np.floor(ls.min() / BIN_WIDTH) * BIN_WIDTH for _, ls in data)
hi = max(np.ceil(ls.max() / BIN_WIDTH) * BIN_WIDTH for _, ls in data)
bins = np.arange(lo, hi + BIN_WIDTH, BIN_WIDTH)

fig, ax = plt.subplots(figsize=(12, 5.6))
out = {"bin_edges_log10": bins.tolist(), "bin_width": BIN_WIDTH, "series": []}
for (label, log_score), color in zip(data, PALETTE):
    counts, edges = np.histogram(log_score, bins=bins, density=True)
    ax.step(edges[:-1], counts, where="post", color=color, linewidth=1.8, label=f"{label} (n={len(log_score)})")
    out["series"].append({"label": label, "n": len(log_score), "density": counts.tolist()})

ax.set_xlabel("log10(IsoRank alignment score)")
ax.set_ylabel("Density (each curve integrates to 1)")
ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
ax.set_axisbelow(True)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.set_title("IsoRank alignment score distributions: paper's rat reference vs. our ncrassa<->calbicans\n"
              "(density-normalized so differing sample sizes are comparable; higher/closer to 0 = stronger match)",
              fontsize=11)
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/rat_vs_ours_overlay.png", dpi=150, bbox_inches="tight")
print("wrote results/rat_vs_ours_overlay.png")

with open("results/rat_vs_ours_overlay.json", "w") as f:
    json.dump(out, f, indent=2)
print("wrote results/rat_vs_ours_overlay.json")
