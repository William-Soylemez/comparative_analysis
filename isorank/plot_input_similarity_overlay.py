#!/usr/bin/env python3
"""
Overlay comparison of the PRE-IsoRank sequence-similarity scores themselves
(the "E" input to compute_pairs/compute_isorank, i.e. the net/<a>-<b>.tsv
files) -- not the post-IsoRank alignment scores (see
plot_rat_vs_ours_overlay.py for that). This is the actual input signal we
control, so it's the fairest comparison of "our metric" vs. "the paper's
metric" before any topology propagation reshapes it.

All three are already min-max normalized to [0, 1] within their own file
(ours explicitly in build_rblast.py / build_rblast_bitscore.py; the paper's
rat-bakers.tsv appears to be normalized the same way, though its repo has
no visible preprocessing script to confirm the exact recipe -- see prior
investigation notes). Linear x-axis (not log) since all three already live
in a bounded [0, 1] range, unlike the post-IsoRank scores which spanned many
orders of magnitude. Density (not raw count) on y so differing pair counts
(rat-bakers: 8606, ours: 30155 each) are comparable.

Usage:
    python3 plot_input_similarity_overlay.py
Expects net/ncrassa-calbicans.tsv, net/ncrassa_bs-calbicans_bs.tsv,
        /tmp/netalign_repo_check/data/intact/rat-bakers.tsv
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2a78d6", "#1baf7a", "#eda100"]
BIN_WIDTH = 0.02
RAT_PAIRS = "/tmp/netalign_repo_check/data/intact/rat-bakers.tsv"

series = [
    ("rat-bakers (paper input)", RAT_PAIRS),
    ("ncrassa-calbicans, E-score input", "net/ncrassa-calbicans.tsv"),
    ("ncrassa-calbicans, bitscore input", "net/ncrassa_bs-calbicans_bs.tsv"),
]

data = []
for label, path in series:
    scores = pd.read_csv(path, sep="\t")["score"].to_numpy()
    data.append((label, scores))

bins = np.arange(0, 1 + BIN_WIDTH, BIN_WIDTH)

fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
out = {"bin_edges": bins.tolist(), "bin_width": BIN_WIDTH, "series": []}
for (label, scores), color in zip(data, PALETTE):
    counts, edges = np.histogram(scores, bins=bins, density=True)
    for ax in axes:
        ax.step(edges[:-1], counts, where="post", color=color, linewidth=1.8,
                 label=f"{label} (n={len(scores)})")
    out["series"].append({"label": label, "n": len(scores), "density": counts.tolist(),
                           "min": float(scores.min()), "max": float(scores.max()),
                           "median": float(np.median(scores))})

axes[0].set_title("Full range [0, 1]")
axes[1].set_yscale("log")
axes[1].set_title("Same data, log y (reveals the tail)")
for ax in axes:
    ax.set_xlabel("Pre-IsoRank sequence-similarity score (min-max normalized within its own file)")
    ax.set_ylabel("Density (each curve integrates to 1)")
    ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
axes[0].legend(frameon=False, fontsize=8, loc="upper right")

fig.suptitle("Sequence-similarity score distributions (pre-IsoRank input): paper's rat-bakers pairs vs. our ncrassa<->calbicans",
             fontsize=11)
fig.tight_layout()
fig.savefig("input_similarity_overlay.png", dpi=150, bbox_inches="tight")
print("wrote input_similarity_overlay.png")

with open("input_similarity_overlay.json", "w") as f:
    json.dump(out, f, indent=2)
print("wrote input_similarity_overlay.json")
