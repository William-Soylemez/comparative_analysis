#!/usr/bin/env python3
"""
Analysis A — Network topology fingerprint.

For each species' predicted PPI network (`*_network.positive.tsv`), compute basic
size-fair topology metrics and plot overlaid degree distributions so the shapes
can be eyeballed side by side.

Usage:
    python3 analysis_a_topology.py [SPECIES_DIR ...]

With no arguments it auto-discovers every species subdirectory that contains a
`*_network.positive.tsv` file. Outputs `topology.json` and
`topology_degree_distribution.png` in the current directory.
"""

import sys
import os
import json
import glob
from collections import Counter

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Categorical palette (dataviz skill, fixed slot order) — one color per species.
PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]


def discover_species(args):
    """Return list of (name, network_path). Args are dirs; empty = auto-discover."""
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_network.positive.tsv")
    )
    out = []
    for d in dirs:
        hits = glob.glob(f"{d.rstrip('/')}/*_network.positive.tsv")
        if not hits:
            print(f"  ! no network file in {d}, skipping", file=sys.stderr)
            continue
        name = os.path.basename(hits[0]).replace("_network.positive.tsv", "")
        out.append((name, hits[0]))
    return out


def load_degrees(path):
    """Read an undirected edge list; return degree Counter and unique edge set."""
    deg = Counter()
    seen = set()
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            if a == b:
                continue  # skip self-loops
            key = (a, b) if a < b else (b, a)
            if key in seen:
                continue
            seen.add(key)
            deg[a] += 1
            deg[b] += 1
    return deg, seen


def metrics(name, path):
    deg, edges = load_degrees(path)
    n_edges = len(edges)
    degs = np.array(list(deg.values()))
    n_nodes = len(deg)
    density = (2 * n_edges) / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0.0
    top = sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:5]

    g = nx.Graph()
    g.add_edges_from(edges)
    global_clustering = nx.transitivity(g)  # 3*triangles / connected triples

    return {
        "species": name,
        "nodes": n_nodes,
        "edges": n_edges,
        "density": density,
        "mean_degree": float(degs.mean()),
        "median_degree": float(np.median(degs)),
        "max_degree": int(degs.max()),
        "global_clustering_coefficient": global_clustering,
        "top_hubs": [{"protein": p, "degree": d} for p, d in top],
        "_degrees": degs,  # kept for plotting, stripped before JSON
    }


def plot(results, outfile):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    for i, r in enumerate(results):
        color = PALETTE[i % len(PALETTE)]
        degs = r["_degrees"]

        # Panel 1: degree histogram (linear x, log y) — the bulk of low-degree nodes.
        ax1.hist(degs, bins=range(1, degs.max() + 2), histtype="step",
                 linewidth=2, color=color, label=r["species"])

        # Panel 2: complementary CDF on log-log — size-fair, reveals scale-free tail.
        vals = np.sort(degs)
        ccdf = 1.0 - np.arange(len(vals)) / len(vals)  # fraction of nodes with degree >= k
        ax2.plot(vals, ccdf, linewidth=2, color=color, label=r["species"])

    ax1.set_yscale("log")
    ax1.set_xlabel("Degree (number of interaction partners)")
    ax1.set_ylabel("Number of proteins (log)")
    ax1.set_title("Degree histogram")
    ax1.grid(True, which="both", color="#e6e6e3", linewidth=0.6)
    ax1.set_axisbelow(True)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Degree k (log)")
    ax2.set_ylabel("Fraction of proteins with degree ≥ k")
    ax2.set_title("Degree distribution (CCDF, size-fair)")
    ax2.grid(True, which="both", color="#e6e6e3", linewidth=0.6)
    ax2.set_axisbelow(True)

    ax2.legend(frameon=False, title="Species")
    fig.suptitle("Analysis A — PPI network topology", fontsize=14, y=0.99)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def plot_bars(results, outfile):
    """Small panels (own scale each) for nodes, edges, density, clustering coefficient."""
    names = [r["species"] for r in results]
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(results))]
    specs = [
        ("nodes", "Proteins in network", "{:,.0f}"),
        ("edges", "Predicted interactions", "{:,.0f}"),
        ("density", "Density (edges / possible)", "{:.4f}"),
        ("global_clustering_coefficient", "Global clustering coeff.", "{:.4f}"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.6))
    x = range(len(names))
    for ax, (key, title, fmt) in zip(axes, specs):
        vals = [r[key] for r in results]
        bars = ax.bar(x, vals, color=colors, width=0.65)
        ax.set_title(title)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=20, ha="right", fontsize=9)
        ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.margins(y=0.15)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, fmt.format(v),
                    ha="center", va="bottom", fontsize=9, color="#52514e")
    fig.suptitle("Analysis A — network size & density", fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No species networks found.", file=sys.stderr)
        sys.exit(1)

    results = []
    print(f"{'species':<20}{'nodes':>8}{'edges':>10}{'density':>10}{'mean_deg':>10}{'max_deg':>9}{'clust_coef':>12}")
    for name, path in species:
        r = metrics(name, path)
        results.append(r)
        print(f"{r['species']:<20}{r['nodes']:>8}{r['edges']:>10}{r['density']:>10.5f}"
              f"{r['mean_degree']:>10.2f}{r['max_degree']:>9}{r['global_clustering_coefficient']:>12.5f}")

    plot(results, "topology_degree_distribution.png")
    plot_bars(results, "topology_size_bars.png")

    # Strip the numpy arrays before serializing.
    for r in results:
        r.pop("_degrees", None)
    with open("topology.json", "w") as f:
        json.dump(results, f, indent=2)
    print("  wrote topology.json")


if __name__ == "__main__":
    main()
