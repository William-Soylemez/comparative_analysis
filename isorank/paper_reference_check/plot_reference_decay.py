#!/usr/bin/env python3
"""Same decay-curve view as isorank/plot_score_histogram.py, applied to the
bakers/rat reference alignment (real paper example data). Uses log10(score),
not -log10(score): score is a "bigger = more confident" similarity, so
log10(score) keeps that direction (closer to 0 = strong match, more
negative = near noise floor) instead of inverting it."""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HUE = "#2a78d6"
BIN_WIDTH = 0.1

df = pd.read_csv("bakers_rat_alignment_scored.tsv", sep="\t")
scores = df["score"].to_numpy()
log_score = np.log10(scores)
n = len(log_score)

lo = np.floor(log_score.min() / BIN_WIDTH) * BIN_WIDTH
hi = np.ceil(log_score.max() / BIN_WIDTH) * BIN_WIDTH
bins = np.arange(lo, hi + BIN_WIDTH, BIN_WIDTH)
counts, edges = np.histogram(log_score, bins=bins)

fig, ax = plt.subplots(figsize=(11, 5.2))
ax.bar(edges[:-1], counts, width=BIN_WIDTH * 0.92, align="edge", color=HUE)
ax.set_xlabel("log10(IsoRank alignment score)")
ax.set_ylabel(f"Count of aligned pairs (bin width {BIN_WIDTH})")
ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
ax.set_axisbelow(True)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.set_title("IsoRank alignment score decay curve, bakers <-> rat (paper's own reference data)\n"
             f"(n={n} pairs, 100% of the smaller species' nodes matched; "
             f"higher/closer to 0 = strong match, more negative = near noise floor)", fontsize=11)
fig.tight_layout()
fig.savefig("bakers_rat_score_decay.png", dpi=150, bbox_inches="tight")
print("wrote bakers_rat_score_decay.png")

out = {"n_pairs": n, "bin_width": BIN_WIDTH,
       "bin_edges_log10": edges.tolist(), "bin_counts": counts.tolist()}
with open("bakers_rat_score_decay.json", "w") as f:
    json.dump(out, f, indent=2)
print("wrote bakers_rat_score_decay.json")
