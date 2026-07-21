#!/usr/bin/env python3
"""
Analysis B (Option 1) — Functional repertoire comparison.

Purely from `*_GO_map.csv` (no network). For each species, count how many
proteins carry each GO term, normalize to a fraction of the proteome (fair across
different proteome sizes), and compare the most frequent / most divergent terms
across species.

GO term *names* are scraped from the `*_human_readable.txt` files so labels are
readable — no external ontology or dependencies required.

Usage:
    python3 analysis_b_go_frequency.py [SPECIES_DIR ...]

Outputs `go_frequency.json`, `go_frequency_heatmap.png`, and
`go_frequency_differences.png` in the current directory.
"""

import sys
import os
import re
import csv
import json
import glob
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]

TOP_N_HEATMAP = 25      # most frequent terms to show in the heatmap
TOP_N_DIFF = 15         # most cross-species-divergent terms to show
DIFF_MIN_FRACTION = 0.02  # a term must reach this in >=1 species to enter the diff view


def discover_species(args):
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_GO_map.csv")
    )
    out = []
    for d in dirs:
        hits = glob.glob(f"{d.rstrip('/')}/*_GO_map.csv")
        if not hits:
            continue
        name = os.path.basename(hits[0]).replace("_GO_map.csv", "")
        out.append((name, d.rstrip("/"), hits[0]))
    return out


def scrape_go_names(species):
    """Build GO_id -> readable name from all human_readable.txt files."""
    names = {}
    pat = re.compile(r"(GO:\d+)\s*-\s*<([^>]+)>")
    for _, d, _ in species:
        for hr in glob.glob(f"{d}/*_human_readable.txt"):
            with open(hr) as f:
                for m in pat.finditer(f.read()):
                    names.setdefault(m.group(1), m.group(2))
    return names


def count_terms(go_map_path):
    """Return (proteome_size, Counter of proteins carrying each GO term)."""
    counts = Counter()
    nprot = 0
    with open(go_map_path) as f:
        r = csv.DictReader(f)
        for row in r:
            nprot += 1
            # set() so we count *proteins carrying* a term, not repeated mentions
            terms = {g for g in (row.get("GO_list") or "").split(";") if g}
            counts.update(terms)
    return nprot, counts


def label(go_id, names, maxlen=42):
    nm = names.get(go_id)
    if not nm:
        return go_id
    if len(nm) > maxlen:
        nm = nm[: maxlen - 1] + "…"
    return f"{nm}"


def plot_heatmap(order, frac, species_names, names, outfile):
    rows = [label(g, names) for g in order]
    mat = np.array([[frac[s].get(g, 0.0) for s in species_names] for g in order])

    fig, ax = plt.subplots(figsize=(9, 0.42 * len(order) + 1.5))
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0)
    ax.set_xticks(range(len(species_names)))
    ax.set_xticklabels(species_names, rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(rows, fontsize=8)
    # annotate each cell with the percentage
    thresh = mat.max() * 0.6
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, f"{100*v:.0f}%", ha="center", va="center", fontsize=7,
                    color="white" if v > thresh else "#333333")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fraction of proteome carrying term")
    ax.set_title(f"Analysis B — {TOP_N_HEATMAP} most frequent GO terms", fontsize=13)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def plot_differences(frac, species_names, names, outfile):
    # candidate terms: reach DIFF_MIN_FRACTION in at least one species
    all_terms = set()
    for s in species_names:
        all_terms.update(frac[s].keys())
    scored = []
    for g in all_terms:
        vals = [frac[s].get(g, 0.0) for s in species_names]
        if max(vals) >= DIFF_MIN_FRACTION:
            scored.append((max(vals) - min(vals), g, vals))
    scored.sort(reverse=True)
    top = scored[:TOP_N_DIFF][::-1]  # reverse so biggest is at top of barh

    y = np.arange(len(top))
    h = 0.8 / len(species_names)
    fig, ax = plt.subplots(figsize=(10, 0.5 * len(top) + 1.5))
    for j, s in enumerate(species_names):
        vals = [100 * t[2][j] for t in top]
        ax.barh(y + j * h, vals, height=h, color=PALETTE[j % len(PALETTE)], label=s)
    ax.set_yticks(y + h * (len(species_names) - 1) / 2)
    ax.set_yticklabels([label(t[1], names) for t in top], fontsize=8)
    ax.set_xlabel("Percent of proteome carrying term")
    ax.grid(True, axis="x", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, title="Species", fontsize=8)
    ax.set_title(f"Analysis B — {TOP_N_DIFF} most divergent GO terms across species", fontsize=13)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No GO_map files found.", file=sys.stderr)
        sys.exit(1)

    go_names = scrape_go_names(species)
    print(f"scraped {len(go_names)} GO term names from human_readable files")

    species_names = []
    proteome = {}
    frac = {}      # species -> {go_id: fraction}
    raw = {}       # species -> {go_id: n_proteins}
    for name, d, path in species:
        nprot, counts = count_terms(path)
        species_names.append(name)
        proteome[name] = nprot
        raw[name] = counts
        frac[name] = {g: c / nprot for g, c in counts.items()}
        print(f"  {name}: {nprot} proteins, {len(counts)} distinct GO terms")

    # top-N most frequent by mean fraction across species (for the heatmap)
    all_terms = set().union(*[set(frac[s]) for s in species_names])
    ranked = sorted(all_terms, key=lambda g: np.mean([frac[s].get(g, 0) for s in species_names]),
                    reverse=True)
    order = ranked[:TOP_N_HEATMAP]

    plot_heatmap(order, frac, species_names, go_names, "go_frequency_heatmap.png")
    plot_differences(frac, species_names, go_names, "go_frequency_differences.png")

    # JSON: proteome sizes + top 100 terms per species with names
    out = {"proteome_sizes": proteome, "species": {}}
    for s in species_names:
        top100 = sorted(raw[s].items(), key=lambda kv: kv[1], reverse=True)[:100]
        out["species"][s] = [
            {"go_id": g, "name": go_names.get(g, ""), "n_proteins": c,
             "fraction": c / proteome[s]}
            for g, c in top100
        ]
    with open("go_frequency.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote go_frequency.json")


if __name__ == "__main__":
    main()
