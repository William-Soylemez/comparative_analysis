#!/usr/bin/env python3
"""
Analysis C on GO-slim categories (repeat of the "across-species dispersion" view,
but using the ~70 high-level slim categories instead of raw GO terms).

Caveat baked in: slim categories are huge (often >1000 proteins), so their absolute
concentration ("fraction in biggest cluster") is expected to be low for everything —
a category that big cannot fit in one ~29-protein cluster. The interesting question
here is therefore RELATIVE: for the same category, does one species concentrate it
more than another?

Usage:
    python3 analysis_c2_goslim_dispersion.py [SPECIES_DIR ...]

Outputs goslim_dispersion_across.png and goslim_dispersion.json.
"""

import sys
import os
import csv
import json
import glob
import math
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from goslim_util import GoSlim

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
OBO = "go_data/go-basic.obo"
SLIM = "go_data/goslim_generic.obo"
TOP_N_ACROSS = 20


def discover_species(args):
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_clusters.json")
    )
    out = []
    for d in dirs:
        d = d.rstrip("/")
        hits = glob.glob(f"{d}/*_clusters.json")
        if hits:
            out.append((os.path.basename(hits[0]).replace("_clusters.json", ""), d))
    return out


def protein_to_cluster(clusters_path):
    p2c = {}
    for cid, c in json.load(open(clusters_path)).items():
        for m in c["members"]:
            p2c[m] = cid
    return p2c


def score_species(name, d, gs):
    p2c = protein_to_cluster(f"{d}/{name}_clusters.json")
    cat_clusters = defaultdict(Counter)
    for row in csv.DictReader(open(f"{d}/{name}_GO_map.csv")):
        cid = p2c.get(row["prot_id"])
        if cid is None:
            continue
        cats = set()
        for t in (row.get("GO_list") or "").split(";"):
            if t:
                cats |= (gs.map_term(t) - ROOTS)
        for c in cats:
            cat_clusters[c][cid] += 1
    scores = {}
    for c, counts in cat_clusters.items():
        n = sum(counts.values())
        if n < 8:
            continue
        frac = max(counts.values()) / n
        H = -sum((k / n) * math.log(k / n) for k in counts.values())
        scores[c] = {"n_proteins": n, "n_clusters": len(counts),
                     "frac_biggest": frac, "norm_entropy": H / math.log(n) if n > 1 else 0.0}
    return scores


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No clusters.json found.", file=sys.stderr)
        sys.exit(1)
    gs = GoSlim(OBO, SLIM)
    species_names = [s[0] for s in species]

    all_scores = {}
    for name, d in species:
        sc = score_species(name, d, gs)
        all_scores[name] = sc
        fr = [v["frac_biggest"] for v in sc.values()]
        print(f"  {name}: {len(sc)} slim categories scored; "
              f"median concentration {np.median(fr):.3f} "
              f"(raw-term C was ~0.11-0.15 — note how much lower this is)")

    shared = set.intersection(*[set(all_scores[sp]) for sp in species_names])
    scored = []
    for c in shared:
        vals = [all_scores[sp][c]["frac_biggest"] for sp in species_names]
        ns = [all_scores[sp][c]["n_proteins"] for sp in species_names]
        scored.append((max(vals) - min(vals), c, vals, ns))
    scored.sort(reverse=True)
    top = scored[:TOP_N_ACROSS][::-1]
    print(f"  {len(shared)} categories shared across all species; showing top {len(top)} by divergence")

    y = np.arange(len(top))
    h = 0.8 / len(species_names)
    fig, ax = plt.subplots(figsize=(10.5, 0.55 * len(top) + 1.6))
    for j, sp in enumerate(species_names):
        bars = ax.barh(y + j * h, [t[2][j] for t in top], height=h,
                       color=PALETTE[j % len(PALETTE)], label=sp)
        for b, t in zip(bars, top):
            ax.text(b.get_width() + 0.003, b.get_y() + b.get_height() / 2,
                    f"n={t[3][j]}", va="center", fontsize=6, color="#777777")
    ax.set_yticks(y + h * (len(species_names) - 1) / 2)
    ax.set_yticklabels([gs.slim_name(t[1]) for t in top], fontsize=8)
    ax.set_xlabel("Concentration (fraction in biggest cluster)")
    ax.grid(True, axis="x", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, title="Species", fontsize=8, loc="lower right")
    ax.set_title(f"Analysis C on GO-slim — {len(top)} most divergent of {len(shared)} categories\n"
                 "(note: absolute concentration is low — slim categories are too big to be modules)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig("goslim_dispersion_across.png", dpi=150, bbox_inches="tight")
    print("  wrote goslim_dispersion_across.png")

    out = {
        "note": "slim categories are large; absolute concentration is low by construction",
        "median_concentration": {sp: float(np.median([v["frac_biggest"] for v in all_scores[sp].values()]))
                                 for sp in species_names},
        "most_divergent": [
            {"go_id": c, "name": gs.slim_name(c),
             "frac_biggest": dict(zip(species_names, vals)),
             "n_proteins": dict(zip(species_names, ns))}
            for rng, c, vals, ns in reversed(top)
        ],
    }
    with open("goslim_dispersion.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote goslim_dispersion.json")


if __name__ == "__main__":
    main()
