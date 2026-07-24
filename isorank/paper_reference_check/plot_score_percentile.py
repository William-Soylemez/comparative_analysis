#!/usr/bin/env python3
"""
Percentile-band view of IsoRank alignment scores: percentile rank on x
(1-percentile-wide bins, full 0-100 range), mean score in that band on y.
Written as both linear and log y-scale versions. Companion to
plot_score_histogram.py's decay-curve view (percentile-rank x-axis vs.
-log10(score) x-axis look at the same data two different ways — this one
shows where a given percentile of confidence sits on the raw score scale,
the decay-curve view shows the actual count distribution).

Usage:
    python3 plot_score_percentile.py <species_a> <species_b> [pct]
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
PCT = int(sys.argv[3]) if len(sys.argv) > 3 else 100
prefix = f"{a}_{b}"

df = pd.read_csv(f"{prefix}_alignment_scored.tsv", sep="\t")
scores = np.sort(df["score"].to_numpy())
n = len(scores)

cutoff_idx = int(np.ceil(PCT / 100 * n)) if PCT < 100 else n
bottom = scores[:cutoff_idx]

bin_means, bin_mins, bin_maxs = [], [], []
for bnd in range(PCT):
    lo = int(np.floor(bnd / 100 * n))
    hi = int(np.floor((bnd + 1) / 100 * n)) if bnd < PCT - 1 else cutoff_idx
    seg = bottom[lo:hi]
    bin_means.append(seg.mean() if len(seg) else np.nan)
    bin_mins.append(seg.min() if len(seg) else np.nan)
    bin_maxs.append(seg.max() if len(seg) else np.nan)

x = np.arange(PCT)
label_bands = set(range(0, PCT, 10)) | {PCT - 1}


def make_plot(yscale, outfile):
    fig, ax = plt.subplots(figsize=(15, 5.2))
    ax.bar(x, bin_means, width=0.85, color=HUE)
    ax.set_xticks(x[::5])
    ax.set_xticklabels([f"{bnd}" for bnd in x[::5]], fontsize=7, rotation=0)
    ax.set_xlabel(f"Percentile of alignment score (0-{PCT}%, full range)")
    ax.set_ylabel("Mean IsoRank score in percentile band"
                  f"{' (log scale)' if yscale == 'log' else ''}")
    if yscale == "log":
        ax.set_yscale("log")
    ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    for bnd in sorted(label_bands):
        v = bin_means[bnd]
        if not np.isnan(v):
            ax.text(bnd, v, f"{v:.2g}", ha="center", va="bottom", fontsize=6,
                    color="#52514e", rotation=90)

    ax.set_title(f"IsoRank alignment scores, {a} <-> {b} — full 0-{PCT}% by percentile band "
                 f"({yscale} y-axis)\n"
                 f"(n={n} pairs total, 100% of the smaller species' nodes matched)", fontsize=11)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"wrote {outfile}")


make_plot("linear", f"{prefix}_score_percentile{PCT}_linear.png")
make_plot("log", f"{prefix}_score_percentile{PCT}_log.png")

out = {
    "species_a": a,
    "species_b": b,
    "n_total_pairs": n,
    "pct": PCT,
    "bin_edges_percentile": list(range(PCT + 1)),
    "bin_mean_score": [float(v) if not np.isnan(v) else None for v in bin_means],
    "bin_min_score": [float(v) if not np.isnan(v) else None for v in bin_mins],
    "bin_max_score": [float(v) if not np.isnan(v) else None for v in bin_maxs],
}
with open(f"{prefix}_score_percentile{PCT}.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"wrote {prefix}_score_percentile{PCT}.json")
