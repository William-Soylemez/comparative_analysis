#!/usr/bin/env python3
"""
Same decay-curve idea as plot_score_histogram.py, but with the raw
(non-log-transformed) IsoRank score on the x-axis, in fixed-width buckets
(default 1e-6). Requested as a "what does it look like without the log
transform" sanity check.

Fair warning baked into the plot itself: because scores span many orders of
magnitude (most pairs are forced matches with near-zero topology-only
support), the overwhelming majority of pairs will land in the very first
bucket [0, bin_width) — that's not a bug, it's the actual shape of the data
before any log transform hides how extreme the skew is. The y-axis (count)
is log-scaled so the sparse right-hand tail (the real matches) stays visible
at all; without that, a linear y-axis would make the tail invisible next to
the near-zero spike.

Usage:
    python3 plot_score_linear.py <species_a> <species_b> [bin_width]
Expects <species_a>_<species_b>_alignment_scored.tsv (from run_full_alignment.py).
"""

import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HUE = "#2a78d6"


def metric_label(name):
    """Similarity metric used to build net/<name>-<other>.tsv, inferred from
    the species alias convention used across this project's scripts: "_bs"
    = build_rblast_bitscore.py (raw bitscore, max across BLAST directions),
    anything else = build_rblast.py (-log10(evalue), min evalue across
    directions, i.e. "E-score")."""
    return "bitscore" if name.endswith("_bs") else "E-score (-log10 evalue)"


a, b = sys.argv[1], sys.argv[2]
BIN_WIDTH = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
prefix = f"{a}_{b}"
METRIC = metric_label(a)

df = pd.read_csv(f"{prefix}_alignment_scored.tsv", sep="\t")
scores = df["score"].to_numpy()
n = len(scores)

hi = np.ceil(scores.max() / BIN_WIDTH) * BIN_WIDTH
bins = np.arange(0, hi + BIN_WIDTH, BIN_WIDTH)
counts, edges = np.histogram(scores, bins=bins)

first_bin_count = int(counts[0])
first_bin_frac = first_bin_count / n

fig, ax = plt.subplots(figsize=(11, 5.2))
ax.bar(edges[:-1], counts, width=BIN_WIDTH * 0.92, align="edge", color=HUE)
ax.set_yscale("log")
ax.set_xlabel(f"IsoRank alignment score (raw, not log-transformed; bucket width {BIN_WIDTH:g})  [similarity metric: {METRIC}]")
ax.set_ylabel("Count of aligned pairs (log scale)")
ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
ax.set_axisbelow(True)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.set_title(f"IsoRank alignment scores, {a} <-> {b} — raw score, not log-transformed — similarity metric: {METRIC}\n"
             f"(n={n} pairs; {first_bin_count} pairs [{100*first_bin_frac:.1f}%] fall in the "
             f"very first bucket [0, {BIN_WIDTH:g}) — count axis is log-scaled so the tail stays visible)",
             fontsize=11)
fig.tight_layout()
outfile = f"{prefix}_score_linear.png"
fig.savefig(outfile, dpi=150, bbox_inches="tight")
print(f"wrote {outfile}")
print(f"  {first_bin_count}/{n} ({100*first_bin_frac:.1f}%) pairs fall in the first bucket [0, {BIN_WIDTH:g})")

out = {
    "species_a": a, "species_b": b, "similarity_metric": METRIC, "n_pairs": n, "bin_width": BIN_WIDTH,
    "first_bin_count": first_bin_count, "first_bin_fraction": first_bin_frac,
    "bin_edges": edges.tolist(), "bin_counts": counts.tolist(),
}
with open(f"{prefix}_score_linear.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"wrote {prefix}_score_linear.json")
