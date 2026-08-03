#!/usr/bin/env python3
"""
Plot n=4 connected-graphlet counts for 3 representative species (spanning phyla
and network density). Two panels:
  A. absolute counts, log scale  -> magnitude/size differences
  B. composition (% of each species' 4-graphlets) -> structural profile

Color = species (categorical slots 1-3 from the validated dataviz palette).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

# --- representative species: sparse -> mid -> dense, 3 different phyla ---
SPECIES = [
    ("GCF_000146465.1", "E. intestinalis\n(Microsporidia, avg deg 39)", "#2a78d6"),
    ("GCF_000182925.2", "N. crassa\n(Ascomycota, avg deg 113)",       "#1baf7a"),
    ("GCF_000149245.1", "C. neoformans\n(Basidiomycota, avg deg 420)", "#eda100"),
]
# graphlet classes ordered sparse->dense (by edge count)
CLASSES = [("P4_path", "P4\npath"), ("claw_star", "claw\n(3-star)"),
           ("C4_cycle", "C4\ncycle"), ("paw", "paw\n(tri+tail)"),
           ("diamond", "diamond\n(K4−e)"), ("K4_clique", "K4\nclique")]

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, SURFACE = "#e1e0d9", "#fcfcfb"

data = {e["accession"]: e for e in json.load(open(os.path.join(RESULTS, "graphlet4_counts.json")))}

counts = np.array([[data[acc]["graphlets"][k] for k, _ in CLASSES]
                   for acc, _, _ in SPECIES], dtype=float)          # [species, class]
totals = counts.sum(axis=1, keepdims=True)
frac = counts / totals * 100.0

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED,
})

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))
x = np.arange(len(CLASSES))
w = 0.26

for i, (acc, label, color) in enumerate(SPECIES):
    off = (i - 1) * w
    axA.bar(x + off, counts[i], w * 0.92, color=color, label=label.split("\n")[0],
            zorder=3)
    axB.bar(x + off, frac[i], w * 0.92, color=color, label=label, zorder=3)
    # relief rule: direct labels on the composition panel
    for xi, f in zip(x + off, frac[i]):
        axB.text(xi, f + 0.6, f"{f:.0f}", ha="center", va="bottom",
                 fontsize=7.5, color=INK2)

# Panel A: absolute, log
axA.set_yscale("log")
axA.set_title("Absolute count of each 4-graphlet  (log scale)",
              fontsize=11, color=INK, loc="left", pad=10)
axA.set_ylabel("number of induced subgraphs")
axA.set_ylim(1e6, 5e11)

# Panel B: composition
axB.set_title("Composition  (% of that species' connected 4-graphlets)",
              fontsize=11, color=INK, loc="left", pad=10)
axB.set_ylabel("share of species' 4-graphlets  (%)")
axB.set_ylim(0, max(frac.max() + 8, 60))

for ax in (axA, axB):
    ax.set_xticks(x)
    ax.set_xticklabels([c for _, c in CLASSES], fontsize=8.5, color=INK2)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(length=0)

axB.legend(frameon=False, fontsize=8.5, loc="upper right",
           labelcolor=INK2, handlelength=1.1)

fig.suptitle("Connected 4-graphlet profiles across the fungal tree",
             fontsize=13, color=INK, x=0.06, ha="left", weight="bold")
fig.text(0.06, 0.005,
         "Exact counts via ORCA. Classes ordered sparse to dense (P4/claw are trees, "
         "K4 is the 4-clique).",
         fontsize=8, color=MUTED, ha="left")
fig.tight_layout(rect=[0, 0.03, 1, 0.94])
os.makedirs(RESULTS, exist_ok=True)
out = os.path.join(RESULTS, "graphlet4_profiles.png")
fig.savefig(out, dpi=200)
print("wrote", out)

# table view (relief-rule mitigation + record)
print("\n%-16s " % "class" + " ".join(f"{lab.split(chr(10))[0]:>16}" for _, lab, _ in SPECIES))
for j, (k, name) in enumerate(CLASSES):
    print("%-16s " % name.replace("\n", " ") +
          " ".join(f"{int(counts[i][j]):>16,}" for i in range(len(SPECIES))))
print("%-16s " % "TOTAL" +
      " ".join(f"{int(totals[i][0]):>16,}" for i in range(len(SPECIES))))
