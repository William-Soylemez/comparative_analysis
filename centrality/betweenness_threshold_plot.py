#!/usr/bin/env python3
"""
Betweenness-centrality threshold picker for one species.

Computes (exact, unweighted) betweenness centrality for every node in a
species' PPI network, then draws the top-tail "decay curve": the betweenness
value at each top-percentile cutoff, from 0.5% to --max-pct in 0.5% steps.
An elbow in that curve is the natural place to draw the "high-centrality" line
(top 1% / 5% / 10% are marked for reference).

Why unweighted: the *_network.positive.tsv edge weights are D-SCRIPT
*similarity* scores, but igraph/networkx treat edge weights as path *lengths*
(distances). Feeding similarities in raw would invert the geometry. Plain
topological betweenness (every edge length 1) is the correct connector-protein
measure and the right quick-and-dirty default. Values are normalized by
(n-1)(n-2)/2 so they are comparable across species of different sizes.

The load / compute / plot steps are exposed as functions so the batch runner
(run_all_betweenness.py) can reuse them.

Usage:
    python3 betweenness_threshold_plot.py SPECIES_DIR [--max-pct 20] [--step 0.5]

Example:
    python3 betweenness_threshold_plot.py ../input/GCF_000146045.2   # S. cerevisiae
"""

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HUE = "#2a78d6"
MARK = "#d64a2a"


def load_graph(species_dir):
    """Read *_network.positive.tsv into a simplified undirected igraph."""
    hits = glob.glob(f"{species_dir.rstrip('/')}/*_network.positive.tsv")
    if not hits:
        sys.exit(f"no *_network.positive.tsv in {species_dir}")
    edges = []
    with open(hits[0]) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != p[1]:
                edges.append((p[0], p[1]))
    g = ig.Graph.TupleList(edges, directed=False)
    g.simplify(multiple=True, loops=True)
    return g


def norm_betweenness(g):
    """Exact, unweighted, size-normalized betweenness as a numpy array."""
    bc = np.array(g.betweenness(directed=False))
    n = g.vcount()
    norm = (n - 1) * (n - 2) / 2.0
    return bc / norm if norm > 0 else bc


def decay_summary(bc, n, max_pct=20.0, step=0.5, marks=(1, 5, 10)):
    """Threshold value at each top-percentile cutoff; returns (pcts, thresholds, dict)."""
    bc_sorted = np.sort(bc)[::-1]
    pcts = np.arange(step, max_pct + step / 2, step)
    ranks = np.clip(np.ceil(pcts / 100.0 * n).astype(int) - 1, 0, n - 1)
    thresholds = bc_sorted[ranks]
    cutoffs = {
        str(p): {
            "n_nodes": int(np.ceil(p / 100.0 * n)),
            "threshold": float(bc_sorted[min(int(np.ceil(p / 100.0 * n)) - 1, n - 1)]),
        } for p in marks if p <= max_pct
    }
    return pcts, thresholds, cutoffs


def plot_decay(acc, g_n, g_m, bc, outdir="results", max_pct=20.0, step=0.5):
    """Write <acc>_betweenness_decay.{png,json}; return the summary dict."""
    os.makedirs(outdir, exist_ok=True)
    n = g_n
    pcts, thresholds, cutoffs = decay_summary(bc, n, max_pct, step)
    bc_sorted = np.sort(bc)[::-1]

    fig, ax = plt.subplots(figsize=(11, 5.4))
    ax.plot(pcts, thresholds, "-o", color=HUE, ms=4, lw=1.6)
    for p in (1, 5, 10):
        if p <= max_pct:
            rank = min(int(np.ceil(p / 100.0 * n)) - 1, n - 1)
            yv = bc_sorted[rank]
            ax.axvline(p, color=MARK, lw=0.9, ls="--", alpha=0.7)
            ax.annotate(f"top {p}%\n{yv:.2e}\n({rank+1} nodes)",
                        xy=(p, yv), xytext=(p + 0.3, yv),
                        fontsize=8, color=MARK, va="bottom")
    ax.set_xlabel("Top percentile of nodes (by betweenness)")
    ax.set_ylabel("Betweenness centrality at cutoff (normalized)")
    ax.grid(True, color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        f"Betweenness threshold decay — {acc}\n"
        f"(n={n:,} nodes, {g_m:,} edges; unweighted; elbow = natural cutoff for "
        f"'high-centrality' set)", fontsize=11)
    fig.tight_layout()
    outpng = os.path.join(outdir, f"{acc}_betweenness_decay.png")
    fig.savefig(outpng, dpi=150, bbox_inches="tight")
    plt.close(fig)

    out = {
        "species": acc,
        "n_nodes": int(n),
        "n_edges": int(g_m),
        "weighted": False,
        "normalized": True,
        "percentiles": pcts.tolist(),
        "thresholds": thresholds.tolist(),
        "top_pct_cutoffs": cutoffs,
    }
    with open(os.path.join(outdir, f"{acc}_betweenness_decay.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("species_dir")
    ap.add_argument("--max-pct", type=float, default=20.0,
                    help="right edge of the percentile sweep (default 20)")
    ap.add_argument("--step", type=float, default=0.5,
                    help="percentile bucket width (default 0.5)")
    ap.add_argument("--outdir", default="results", help="where to write outputs (default results/)")
    args = ap.parse_args()

    acc = os.path.basename(args.species_dir.rstrip("/"))
    t0 = time.time()
    g = load_graph(args.species_dir)
    n, m = g.vcount(), g.ecount()
    print(f"species: {acc}")
    print(f"graph: {n:,} nodes, {m:,} edges  (loaded in {time.time()-t0:.1f}s)")

    t0 = time.time()
    bc = norm_betweenness(g)
    print(f"betweenness computed in {time.time()-t0:.1f}s")
    print(f"betweenness: min={bc.min():.3e} median={np.median(bc):.3e} "
          f"max={bc.max():.3e}")

    out = plot_decay(acc, n, m, bc, outdir=args.outdir, max_pct=args.max_pct, step=args.step)
    print(f"wrote {os.path.join(args.outdir, acc + '_betweenness_decay.png')}")
    print(f"wrote {os.path.join(args.outdir, acc + '_betweenness_decay.json')}")
    c = out["top_pct_cutoffs"]["5"]
    print(f"top 5% cutoff: {c['n_nodes']} nodes, threshold {c['threshold']:.3e}")


if __name__ == "__main__":
    main()
