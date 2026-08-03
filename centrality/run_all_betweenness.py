#!/usr/bin/env python3
"""
Run betweenness-centrality threshold analysis across all 11 local fungal species.

For each species it computes exact, unweighted, size-normalized betweenness
(same method as betweenness_threshold_plot.py), writes the per-species decay
png + json, and records the top-5% cutoff. It then draws one overlay plot of
all 11 decay curves and a combined summary json, so you can confirm the elbow
lands in a consistent place before committing to the top-5% "high-centrality"
definition.

Progress: an overall bar advances one step per finished species; each species
also shows a live elapsed timer while its (blocking, C-level) betweenness call
runs in a worker process, so a slow species never looks hung.

Usage (from this directory):
    ../isorank/venv/bin/python run_all_betweenness.py
    ../isorank/venv/bin/python run_all_betweenness.py --max-pct 20 --step 0.5

By default it discovers ../input/GC[AF]_* species dirs that contain a network file.
"""

import argparse
import glob
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm

from betweenness_threshold_plot import load_graph, norm_betweenness, plot_decay, decay_summary

# accession -> short organism label for the overlay legend
NAMES = {
    "GCF_000146045.2": "S. cerevisiae",
    "GCF_000149245.1": "C. neoformans",
    "GCF_000203795.2": "B. dendrobatidis",
    "GCF_000146465.1": "E. intestinalis",
    "GCF_026210795.1": "R. irregularis",
    "GCF_025024165.1": "K. alabastrina",
    "GCA_014872705.1": "A. bisporus",
    "GCF_000182965.3": "C. albicans",
    "GCF_000182925.2": "N. crassa",
    "GCA_025594325.1": "B. emersonii",
    "GCF_025094135.1": "M. mucedo",
}


def _worker(species_dir):
    """Load graph + compute normalized betweenness in a child process.
    Returns (n, m, bc_list) so the blocking C call can't freeze the timer."""
    g = load_graph(species_dir)
    bc = norm_betweenness(g)
    return g.vcount(), g.ecount(), bc.tolist()


def discover(parent):
    dirs = []
    for d in sorted(glob.glob(os.path.join(parent, "GC[AF]_*"))):
        if os.path.isdir(d) and glob.glob(os.path.join(d, "*_network.positive.tsv")):
            dirs.append(d)
    return dirs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parent", default="../input", help="dir holding the species dirs (default ../input)")
    ap.add_argument("--max-pct", type=float, default=20.0)
    ap.add_argument("--step", type=float, default=0.5)
    ap.add_argument("--outdir", default="results", help="where to write outputs (default results/)")
    args = ap.parse_args()

    species_dirs = discover(args.parent)
    if not species_dirs:
        raise SystemExit(f"no species dirs with a network file under {args.parent}")
    os.makedirs(args.outdir, exist_ok=True)

    results = {}
    grand_t0 = time.time()

    outer = tqdm(species_dirs, desc="species", position=0, unit="sp")
    for sdir in outer:
        acc = os.path.basename(sdir.rstrip("/"))
        label = NAMES.get(acc, acc)
        outer.set_postfix_str(label)

        t0 = time.time()
        # one-shot worker process so the blocking betweenness call can't stall the timer
        with ProcessPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(_worker, sdir)
            with tqdm(total=0, position=1, leave=False,
                      bar_format="  {desc} [{elapsed}]") as timer:
                timer.set_description(f"{label}: load + betweenness")
                while not fut.done():
                    time.sleep(0.25)
                    timer.refresh()
            n, m, bc_list = fut.result()
        bc = np.array(bc_list)
        dt = time.time() - t0

        out = plot_decay(acc, n, m, bc, outdir=args.outdir,
                         max_pct=args.max_pct, step=args.step)
        results[acc] = {**out, "label": label, "seconds": round(dt, 1)}
        c5 = out["top_pct_cutoffs"]["5"]
        tqdm.write(f"  ✓ {label:<18} {n:>6,} nodes  {m:>10,} edges  "
                   f"{dt:5.1f}s   top5%={c5['n_nodes']} nodes @ {c5['threshold']:.2e}")
    outer.close()

    # ---- overlay of all decay curves ----
    fig, ax = plt.subplots(figsize=(12, 6.4))
    cmap = plt.get_cmap("tab20")
    for i, acc in enumerate(results):
        r = results[acc]
        ax.plot(r["percentiles"], r["thresholds"], "-", lw=1.7, color=cmap(i % 20),
                label=f"{r['label']} (n={r['n_nodes']:,})")
    ax.axvline(5, color="#555", lw=1.0, ls="--", alpha=0.7)
    ax.text(5.1, ax.get_ylim()[1] * 0.96, "top 5%", fontsize=9, color="#555", va="top")
    ax.set_xlabel("Top percentile of nodes (by betweenness)")
    ax.set_ylabel("Betweenness centrality at cutoff (normalized)")
    ax.grid(True, color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Betweenness threshold decay across 11 fungal species\n"
                 "(exact, unweighted, size-normalized; dashed line = chosen top-5% cutoff)",
                 fontsize=12)
    ax.legend(fontsize=8, ncol=2, frameon=False)
    fig.tight_layout()
    overlay = os.path.join(args.outdir, "all_species_betweenness_decay.png")
    fig.savefig(overlay, dpi=150, bbox_inches="tight")
    plt.close(fig)

    summ = {
        acc: {
            "label": r["label"], "n_nodes": r["n_nodes"], "n_edges": r["n_edges"],
            "seconds": r["seconds"], "top_pct_cutoffs": r["top_pct_cutoffs"],
        } for acc, r in results.items()
    }
    with open(os.path.join(args.outdir, "all_species_betweenness_summary.json"), "w") as f:
        json.dump(summ, f, indent=2)

    print(f"\ndone: {len(results)} species in {time.time()-grand_t0:.1f}s")
    print(f"wrote {overlay}")
    print(f"wrote {os.path.join(args.outdir, 'all_species_betweenness_summary.json')}")
    print("per-species: <acc>_betweenness_decay.png / .json")


if __name__ == "__main__":
    main()
