#!/usr/bin/env python3
"""
Analysis B, Option 2 — Functional repertoire at GO-slim level.

Same idea as analysis_b_go_frequency.py, but every specific GO term is first rolled
up to high-level GO-slim categories (goslim_generic, pinned to the 2026-03-25 GO
release matching the pipeline). This collapses the ~10k redundant/artifact-laden
raw terms into ~70 interpretable categories, so "how much of each function does each
species have" finally means something.

A protein counts toward a category if any of its GO terms rolls up to it (via
is_a + part_of). Counts are normalized to fraction of the proteome.

Usage:
    python3 analysis_b2_goslim.py [SPECIES_DIR ...]

Outputs goslim_frequency.json, goslim_frequency_heatmap.png,
goslim_frequency_differences.png.
"""

import sys
import os
import csv
import json
import glob
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from goslim_util import GoSlim

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}  # the 3 aspect roots — uninformative
OBO = "go_data/go-basic.obo"
SLIM = "go_data/goslim_generic.obo"
TOP_N_HEATMAP = 20
TOP_N_DIFF = 15


def discover_species(args):
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_GO_map.csv")
    )
    out = []
    for d in dirs:
        d = d.rstrip("/")
        hits = glob.glob(f"{d}/*_GO_map.csv")
        if hits:
            out.append((os.path.basename(hits[0]).replace("_GO_map.csv", ""), hits[0]))
    return out


def count_categories(go_map_path, gs):
    """Return (proteome_size, Counter of proteins per slim category, coverage)."""
    counts = Counter()
    nprot = 0
    covered = 0
    for row in csv.DictReader(open(go_map_path)):
        nprot += 1
        cats = set()
        for t in (row.get("GO_list") or "").split(";"):
            if t:
                cats |= (gs.map_term(t) - ROOTS)
        if cats:
            covered += 1
        counts.update(cats)
    return nprot, counts, covered


def plot_heatmap(order, frac, species_names, gs, outfile):
    rows = [gs.slim_name(c) for c in order]
    mat = np.array([[frac[s].get(c, 0.0) for s in species_names] for c in order])
    fig, ax = plt.subplots(figsize=(9, 0.44 * len(order) + 1.5))
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0)
    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names, rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(rows, fontsize=8)
    thresh = mat.max() * 0.6
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{100*mat[i,j]:.0f}%", ha="center", va="center", fontsize=7,
                    color="white" if mat[i, j] > thresh else "#333333")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fraction of proteome in category")
    ax.set_title(f"Analysis B (GO-slim) — {TOP_N_HEATMAP} most common functional categories",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def plot_differences(frac, species_names, gs, outfile):
    all_cats = set().union(*[set(frac[s]) for s in species_names])
    scored = []
    for c in all_cats:
        vals = [frac[s].get(c, 0.0) for s in species_names]
        if max(vals) >= 0.03:
            scored.append((max(vals) - min(vals), c, vals))
    scored.sort(reverse=True)
    top = scored[:TOP_N_DIFF][::-1]
    y = np.arange(len(top))
    h = 0.8 / len(species_names)
    fig, ax = plt.subplots(figsize=(10, 0.5 * len(top) + 1.6))
    for j, s in enumerate(species_names):
        ax.barh(y + j * h, [100 * t[2][j] for t in top], height=h,
                color=PALETTE[j % len(PALETTE)], label=s)
    ax.set_yticks(y + h * (len(species_names) - 1) / 2)
    ax.set_yticklabels([gs.slim_name(t[1]) for t in top], fontsize=8)
    ax.set_xlabel("Percent of proteome in category")
    ax.grid(True, axis="x", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, title="Species", fontsize=8)
    ax.set_title(f"Analysis B (GO-slim) — {TOP_N_DIFF} most divergent categories", fontsize=13)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No GO_map files found.", file=sys.stderr)
        sys.exit(1)
    gs = GoSlim(OBO, SLIM)
    print(f"loaded {len(gs.slim)} slim categories from {SLIM}")

    species_names, proteome, frac, raw = [], {}, {}, {}
    for name, path in species:
        nprot, counts, covered = count_categories(path, gs)
        species_names.append(name)
        proteome[name] = nprot
        raw[name] = counts
        frac[name] = {c: n / nprot for c, n in counts.items()}
        print(f"  {name}: {nprot} proteins, {covered} ({100*covered/nprot:.0f}%) mapped, "
              f"{len(counts)} categories")

    all_cats = set().union(*[set(frac[s]) for s in species_names])
    ranked = sorted(all_cats,
                    key=lambda c: np.mean([frac[s].get(c, 0) for s in species_names]),
                    reverse=True)
    order = ranked[:TOP_N_HEATMAP]

    plot_heatmap(order, frac, species_names, gs, "goslim_frequency_heatmap.png")
    plot_differences(frac, species_names, gs, "goslim_frequency_differences.png")

    out = {
        "slim": SLIM, "go_release": "2026-03-25",
        "proteome_sizes": proteome,
        "species": {
            s: [
                {"go_id": c, "name": gs.slim_name(c), "aspect": gs.namespace.get(c, ""),
                 "n_proteins": raw[s][c], "fraction": raw[s][c] / proteome[s]}
                for c, _ in raw[s].most_common()
            ]
            for s in species_names
        },
    }
    with open("goslim_frequency.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote goslim_frequency.json")


if __name__ == "__main__":
    main()
