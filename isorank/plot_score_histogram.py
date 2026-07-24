#!/usr/bin/env python3
"""
Decay-curve histogram of IsoRank alignment scores: log10(score) on x (fixed
0.1-wide bins), raw pair COUNT on y — not percentile-rank on x / mean score
on y (the earlier version of this script). Averaging scores within percentile
bins smoothed away the actual shape of the distribution; a plain count
histogram over the natural log scale shows it directly.

Uses log10(score), NOT -log10(score): `score` is already a "bigger = more
confident match" similarity value, not a p-value, so log10(score) preserves
that direction (higher/closer-to-zero = more confident, more negative =
weak/near-noise-floor) instead of inverting it the way -log10 would.

Expect most pairs to sit at very negative log10(score) (near the noise
floor, since IsoRank forces every node of the smaller network to get some
match) with a short rising tail approaching 0 (a handful of strong matches)
— the left edge of that tail is where a confidence cutoff should be chosen.

Usage:
    python3 plot_score_histogram.py <species_a> <species_b> [bin_width]
Expects <species_a>_<species_b>_alignment_scored.tsv (from run_full_alignment.py).
"""

import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HUE = "#2a78d6"  # sequential single-hue, magnitude encoding — one series, no legend needed

a, b = sys.argv[1], sys.argv[2]
BIN_WIDTH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
prefix = f"{a}_{b}"

df = pd.read_csv(f"{prefix}_alignment_scored.tsv", sep="\t")
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
ax.set_title(f"IsoRank alignment score decay curve, {a} <-> {b}\n"
             f"(n={n} pairs, 100% of the smaller species' nodes matched; "
             f"higher/closer to 0 = strong match, more negative = near noise floor)",
             fontsize=11)
fig.tight_layout()
outfile = f"{prefix}_score_decay.png"
fig.savefig(outfile, dpi=150, bbox_inches="tight")
print(f"wrote {outfile}")

out = {
    "species_a": a,
    "species_b": b,
    "n_pairs": n,
    "bin_width": BIN_WIDTH,
    "bin_edges_log10": edges.tolist(),
    "bin_counts": counts.tolist(),
}
with open(f"{prefix}_score_decay.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"wrote {prefix}_score_decay.json")
