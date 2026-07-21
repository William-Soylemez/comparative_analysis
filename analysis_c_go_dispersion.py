#!/usr/bin/env python3
"""
Analysis C — GO modularity: co-localized vs. dispersed.

For each GO term in each species, measure how concentrated its proteins are within
the cluster structure:
  - frac_biggest : fraction of the term's clustered proteins that fall in its single
                   most common cluster (0 = fully dispersed, 1 = one tidy module).
  - norm_entropy : normalized Shannon entropy over occupied clusters
                   (0 = one cluster, 1 = evenly scattered) — a second opinion.

Two views (built from the same per-(species, term) score):
  Q1 (within-species) : distribution of frac_biggest across all terms, per species
                        -> is one species' functional organization more modular?
  Q2 (across-species) : same term compared across species
                        -> which specific functions group in some species but not others?

Usage:
    python3 analysis_c_go_dispersion.py [SPECIES_DIR ...]

Outputs go_dispersion.json, go_dispersion_within.png, go_dispersion_across.png.
"""

import sys
import os
import re
import csv
import json
import glob
import math
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]

MIN_PROTEINS = 8     # a term needs >= this many clustered proteins to be scored (Q1)
MIN_ACROSS = 15      # for Q2, require >= this many proteins in EVERY species (kills threshold noise)
TOP_N_ACROSS = 15    # most cross-species-divergent terms to show in Q2


def discover_species(args):
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_clusters.json")
    )
    out = []
    for d in dirs:
        d = d.rstrip("/")
        hits = glob.glob(f"{d}/*_clusters.json")
        if not hits:
            continue
        name = os.path.basename(hits[0]).replace("_clusters.json", "")
        out.append((name, d))
    return out


def scrape_go_names(species):
    names = {}
    pat = re.compile(r"(GO:\d+)\s*-\s*<([^>]+)>")
    for _, d in species:
        for hr in glob.glob(f"{d}/*_human_readable.txt"):
            with open(hr) as f:
                for m in pat.finditer(f.read()):
                    names.setdefault(m.group(1), m.group(2))
    return names


def protein_to_cluster(clusters_path):
    d = json.load(open(clusters_path))
    p2c = {}
    for cid, c in d.items():
        for m in c["members"]:
            p2c[m] = cid
    return p2c


def score_species(name, d):
    """Return {go_id: dict(scores)} for terms with >= MIN_PROTEINS clustered proteins."""
    p2c = protein_to_cluster(f"{d}/{name}_clusters.json")
    # term -> Counter over cluster ids (only proteins that are clustered)
    term_clusters = defaultdict(Counter)
    with open(f"{d}/{name}_GO_map.csv") as f:
        for row in csv.DictReader(f):
            pid = row["prot_id"]
            cid = p2c.get(pid)
            if cid is None:
                continue  # protein not in any cluster
            terms = {g for g in (row.get("GO_list") or "").split(";") if g}
            for g in terms:
                term_clusters[g][cid] += 1

    scores = {}
    for g, counts in term_clusters.items():
        n = sum(counts.values())
        if n < MIN_PROTEINS:
            continue
        frac_biggest = max(counts.values()) / n
        # normalized Shannon entropy over occupied clusters
        H = -sum((c / n) * math.log(c / n) for c in counts.values())
        norm_entropy = H / math.log(n) if n > 1 else 0.0
        scores[g] = {
            "n_proteins": n,
            "n_clusters": len(counts),
            "frac_biggest": frac_biggest,
            "norm_entropy": norm_entropy,
        }
    return scores


def plot_within(all_scores, species_names, outfile):
    """Q1: distribution of frac_biggest per species (violin)."""
    data = [[s["frac_biggest"] for s in all_scores[sp].values()] for sp in species_names]
    fig, ax = plt.subplots(figsize=(8, 5.2))
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(PALETTE[i % len(PALETTE)])
        body.set_edgecolor(PALETTE[i % len(PALETTE)])
        body.set_alpha(0.55)
    parts["cmedians"].set_color("#333333")
    for i, vals in enumerate(data):
        ax.text(i + 1, np.median(vals), f"  med {np.median(vals):.2f}",
                va="center", fontsize=9, color="#333333")
    ax.set_xticks(range(1, len(species_names) + 1))
    ax.set_xticklabels(species_names, rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("Concentration (fraction in biggest cluster)")
    ax.set_ylim(0, 1.02)
    ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Analysis C, Q1 — how modular is each species?\n"
                 "(higher = GO terms concentrate into single clusters)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def plot_across(all_scores, species_names, names, outfile):
    """Q2: same term across species — most divergent by frac_biggest.

    Robustness filters: term must have >= MIN_ACROSS proteins in EVERY species
    (removes coarse threshold noise) and not be an 'obsolete' GO term.
    """
    shared = set.intersection(*[set(all_scores[sp]) for sp in species_names])
    scored = []
    for g in shared:
        name = names.get(g, "")
        if all(all_scores[sp][g]["n_proteins"] >= MIN_ACROSS for sp in species_names) \
           and name and not name.lower().startswith("obsolete"):
            vals = [all_scores[sp][g]["frac_biggest"] for sp in species_names]
            ns = [all_scores[sp][g]["n_proteins"] for sp in species_names]
            scored.append((max(vals) - min(vals), g, vals, ns))
    scored.sort(reverse=True)
    print(f"  Q2 funnel: {len(shared)} shared (>= {MIN_PROTEINS} in all) "
          f"-> {len(scored)} candidates (>= {MIN_ACROSS} in all, named, non-obsolete)")
    # Collapse GO terms that are perfectly correlated across species (identical
    # protein counts + concentration in every species) — these ride the same
    # Pfam domain / protein set, so they are one signal, not many.
    seen = {}
    dedup = []
    for rng, g, vals, ns in scored:
        sig = (tuple(ns), tuple(round(v, 3) for v in vals))
        if sig in seen:
            seen[sig] += 1
            continue
        seen[sig] = 1
        dedup.append((rng, g, vals, ns))
    collapsed = len(scored) - len(dedup)
    n_distinct = len(dedup)
    top = dedup[:TOP_N_ACROSS][::-1]
    print(f"  Q2: {n_distinct} distinct candidates after collapsing "
          f"{collapsed} redundant (same-domain) terms")

    y = np.arange(len(top))
    h = 0.8 / len(species_names)
    fig, ax = plt.subplots(figsize=(10.5, 0.55 * len(top) + 1.6))
    for j, sp in enumerate(species_names):
        vals = [t[2][j] for t in top]
        bars = ax.barh(y + j * h, vals, height=h, color=PALETTE[j % len(PALETTE)], label=sp)
        for b, t in zip(bars, top):
            ax.text(b.get_width() + 0.01, b.get_y() + b.get_height() / 2,
                    f"n={t[3][j]}", va="center", fontsize=6, color="#777777")
    labels = [(names.get(t[1], t[1])[:44]) for t in top]
    ax.set_yticks(y + h * (len(species_names) - 1) / 2)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Concentration (fraction in biggest cluster)")
    ax.set_xlim(0, 1.15)
    ax.grid(True, axis="x", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, title="Species", fontsize=8, loc="lower right")
    ax.set_title(f"Analysis C, Q2 — {len(top)} most divergent of {n_distinct} candidate terms\n"
                 f"(named, non-obsolete, ≥{MIN_ACROSS} proteins in every species; "
                 f"ranked by max−min concentration)", fontsize=11)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")
    return top


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No clusters.json files found.", file=sys.stderr)
        sys.exit(1)
    go_names = scrape_go_names(species)

    species_names = [s[0] for s in species]
    all_scores = {}
    for name, d in species:
        sc = score_species(name, d)
        all_scores[name] = sc
        fr = [v["frac_biggest"] for v in sc.values()]
        print(f"  {name}: {len(sc)} terms scored (>={MIN_PROTEINS} proteins); "
              f"median concentration {np.median(fr):.2f}")

    plot_within(all_scores, species_names, "go_dispersion_within.png")
    top_across = plot_across(all_scores, species_names, go_names, "go_dispersion_across.png")

    # JSON output
    out = {
        "params": {"min_proteins": MIN_PROTEINS},
        "within_species": {
            sp: {
                "median_frac_biggest": float(np.median([v["frac_biggest"] for v in all_scores[sp].values()])),
                "n_terms_scored": len(all_scores[sp]),
            }
            for sp in species_names
        },
        "per_species_terms": {
            sp: {g: {**v, "name": go_names.get(g, "")} for g, v in all_scores[sp].items()}
            for sp in species_names
        },
        "most_divergent_terms": [
            {"go_id": g, "name": go_names.get(g, ""), "range": rng,
             "frac_biggest": dict(zip(species_names, vals)),
             "n_proteins": dict(zip(species_names, ns))}
            for rng, g, vals, ns in reversed(top_across)
        ],
    }
    with open("go_dispersion.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote go_dispersion.json")


if __name__ == "__main__":
    main()
