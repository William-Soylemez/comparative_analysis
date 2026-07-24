#!/usr/bin/env python3
"""
Heatmap of GO enrichment in high-betweenness ("connector") proteins across the
11 fungal species. Rows = GO categories, columns = species (grouped by phylum
from early-diverging on the right to Dikarya on the left), cell = fold
enrichment, shown only where the term is significant (q <= --qmax); other cells
are left blank so the significant pattern is what you read.

Sequential single-hue (blue) magnitude encoding — one measure (fold), so no
legend, just a colorbar; every shown cell is also direct-labeled with its fold
value, so magnitude never rests on color alone.

Usage (from this directory):
    ../isorank/venv/bin/python plot_enrichment_heatmap.py slim
    ../isorank/venv/bin/python plot_enrichment_heatmap.py direct --min-species 4 --max-rows 30
Reads <acc>_enrichment_<tag>.tsv (from fisher_enrichment.py).
"""

import argparse
import csv
import glob
import json
import statistics

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# blue sequential ramp (dataviz reference palette, steps 100 -> 700)
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SURFACE = "#fcfcfb"
BLANK = "#f0efec"  # neutral for non-significant cells

# species columns, grouped Dikarya -> early-diverging, with phylum for separators
COLUMN_ORDER = [
    ("GCF_000146045.2", "S. cerevisiae", "Ascomycota"),
    ("GCF_000182965.3", "C. albicans", "Ascomycota"),
    ("GCF_000182925.2", "N. crassa", "Ascomycota"),
    ("GCF_000149245.1", "C. neoformans", "Basidiomycota"),
    ("GCA_014872705.1", "A. bisporus", "Basidiomycota"),
    ("GCF_026210795.1", "R. irregularis", "Mucoromycota"),
    ("GCF_025094135.1", "M. mucedo", "Mucoromycota"),
    ("GCF_025024165.1", "K. alabastrina", "Zoopago."),
    ("GCF_000203795.2", "B. dendrobatidis", "Chytrid."),
    ("GCA_025594325.1", "B. emersonii", "Blasto."),
    ("GCF_000146465.1", "E. intestinalis", "Micro."),
]


def load(tag, qmax):
    """Return {acc: {name: (fold, q)}} for all species with a <tag> TSV."""
    data = {}
    for tsv in glob.glob(f"*_enrichment_{tag}.tsv"):
        acc = tsv.split("_enrichment")[0]
        d = {}
        with open(tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                d[row["name"]] = (float(row["fold_enrichment"]), float(row["q"]))
        data[acc] = d
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag", help="enrichment tag, e.g. slim, direct, degree_slim, degree_direct")
    ap.add_argument("--qmax", type=float, default=0.05)
    ap.add_argument("--min-species", type=int, default=2,
                    help="keep a category if significant in >= this many species")
    ap.add_argument("--max-rows", type=int, default=28)
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="accession(s) to drop from the columns, e.g. GCF_000146465.1")
    args = ap.parse_args()

    data = load(args.tag, args.qmax)
    cols = [(a, n, p) for (a, n, p) in COLUMN_ORDER
            if a in data and a not in args.exclude]

    # rank categories by # species significant, then median fold (shown cols only)
    sig_folds = {}
    for acc, _, _ in cols:
        for name, (fold, q) in data[acc].items():
            if q <= args.qmax:
                sig_folds.setdefault(name, []).append(fold)
    rows = [(name, len(f), statistics.median(f)) for name, f in sig_folds.items()
            if len(f) >= args.min_species]
    rows.sort(key=lambda r: (-r[1], -r[2]))
    rows = rows[:args.max_rows]
    row_names = [r[0] for r in rows]

    # build matrix of fold (NaN where not significant / absent)
    M = np.full((len(row_names), len(cols)), np.nan)
    for i, name in enumerate(row_names):
        for j, (acc, _, _) in enumerate(cols):
            fold, q = data[acc].get(name, (np.nan, 1.0))
            if q <= args.qmax and not np.isnan(fold):
                M[i, j] = fold

    vmax = float(np.nanpercentile(M, 95)) if np.isfinite(M).any() else 3.0
    vmax = max(2.0, round(vmax, 1))
    cmap = LinearSegmentedColormap.from_list("blues", BLUE)
    cmap.set_bad(BLANK)

    fig_h = 1.6 + 0.34 * len(row_names)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=1.0, vmax=vmax)

    # 2px surface gaps between cells
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_names), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)

    # annotate significant cells with fold; text color by cell luminance
    for i in range(len(row_names)):
        for j in range(len(cols)):
            v = M[i, j]
            if not np.isnan(v):
                frac = min(1.0, (v - 1.0) / (vmax - 1.0))
                tc = "#ffffff" if frac > 0.55 else INK
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                        fontsize=7.5, color=tc)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([n for _, n, _ in cols], rotation=40, ha="right",
                       fontsize=8.5, style="italic", color=INK2)
    ax.set_yticks(range(len(row_names)))
    ax.set_yticklabels(row_names, fontsize=8.5, color=INK)

    # phylum group separators (thicker lines between phyla)
    phyla = [p for _, _, p in cols]
    for j in range(1, len(cols)):
        if phyla[j] != phyla[j - 1]:
            ax.axvline(j - 0.5, color=MUTED, lw=1.4)
    # phylum labels along the top
    j = 0
    while j < len(cols):
        k = j
        while k + 1 < len(cols) and phyla[k + 1] == phyla[j]:
            k += 1
        single = (k == j)
        ax.text((j + k) / 2, -0.7, phyla[j], ha="center", va="bottom",
                fontsize=7.5, color=MUTED, rotation=(30 if single else 0))
        j = k + 1

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, extend="max")
    cb.set_label("fold enrichment", fontsize=8.5, color=INK2)
    cb.ax.tick_params(labelsize=7.5, color=MUTED, labelcolor=INK2)
    cb.outline.set_visible(False)

    tag_label = ("GO-slim categories" if "slim" in args.tag
                 else "propagated GO terms" if "propagated" in args.tag
                 else "raw GO terms")
    centrality = "degree" if args.tag.startswith("degree") else "betweenness"
    n_sig_rows = len(row_names)
    ax.set_title(
        f"Functional enrichment of high-{centrality} connector proteins "
        f"({tag_label})\n"
        f"top 5% of each species' network by {centrality} vs annotated background; "
        f"cells shown where q≤{args.qmax} (Fisher, BH-FDR); blank = not significant",
        fontsize=11, color=INK, pad=26)
    fig.tight_layout()
    out = f"enrichment_{args.tag}_heatmap.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out}  ({n_sig_rows} categories x {len(cols)} species)")


if __name__ == "__main__":
    main()
