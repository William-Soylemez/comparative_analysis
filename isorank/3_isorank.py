#!/usr/bin/env python3
"""
Step 3/3: IsoRank alignment for a pair of species.

Blends network topology (net/<a>.tsv, net/<b>.tsv) with the normalized-bitscore
sequence similarity (net/<a>-<b>.tsv) via IsoRank, then greedily assigns a
one-to-one matching covering 100% of the SMALLER network's nodes, recording the
IsoRank score at every matched pair. No score cutoff is applied here -- keeping
every pair + its score lets a confidence threshold be chosen downstream.

Calls netalign's compute_adjacency / compute_pairs / compute_isorank directly
(rather than the `netalign isorank` CLI) because the CLI's greedy assignment
returns only index pairs, not the per-match score we want to keep.

Usage: python3 3_isorank.py <a> <b>
  e.g. python3 3_isorank.py scer calbicans
Expects: net/<a>.tsv, net/<b>.tsv, net/<a>-<b>.tsv (from 2_similarity.py)
Writes:  results/<a>_<b>_alignment_scored.tsv (one row per node of the smaller
             network: a-protein, b-protein, score)
         results/<a>_<b>_alignment_scores.json (summary stats)
"""

import json
import os
import sys

import numpy as np
import pandas as pd

from netalign.approx_isorank.io_utils import compute_adjacency, compute_pairs
from netalign.approx_isorank.isorank_compute import compute_isorank

ALPHA = 0.7   # IsoRank topology/sequence blend (0 = sequence only, 1 = topology only)
NITER = 20

HERE = os.path.dirname(os.path.abspath(__file__))
NET = os.path.join(HERE, "net")
RESULTS = os.path.join(HERE, "results")


def greedy_match(R):
    """One-to-one greedy assignment over the IsoRank matrix R, covering all of
    the smaller dimension; returns [(row_idx, col_idx, score), ...]."""
    work = R.copy()
    aligned = []
    for _ in range(min(R.shape)):
        maxcols = np.argmax(work, axis=1)
        row_best = work[np.arange(work.shape[0]), maxcols]
        r = int(np.argmax(row_best))
        c = int(maxcols[r])
        aligned.append((r, c, float(work[r, c])))
        work[:, c] = -1
        work[r, :] = -1
    return aligned


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python3 3_isorank.py <a> <b>")
    a, b = sys.argv[1], sys.argv[2]
    os.makedirs(RESULTS, exist_ok=True)

    print(f"Loading networks and reciprocal-BLAST pairs ({a} <-> {b}) ...")
    df1 = pd.read_csv(os.path.join(NET, f"{a}.tsv"), sep="\t", header=None)
    df2 = pd.read_csv(os.path.join(NET, f"{b}.tsv"), sep="\t", header=None)
    dpairs = pd.read_csv(os.path.join(NET, f"{a}-{b}.tsv"), sep="\t")

    Af1, nA1 = compute_adjacency(df1)
    Af2, nA2 = compute_adjacency(df2)
    print(f"  {a}: {len(nA1)} nodes, {b}: {len(nA2)} nodes")

    E = compute_pairs(dpairs, nA1, nA2, a, b)

    print("Computing IsoRank similarity matrix ...")
    R = compute_isorank(Af1, Af2, E, alpha=ALPHA, maxiter=NITER)[-1]
    print(f"  R shape: {R.shape}")

    print("Greedy one-to-one assignment (100% of the smaller network) ...")
    aligned = greedy_match(R)
    rnA1 = {v: k for k, v in nA1.items()}
    rnA2 = {v: k for k, v in nA2.items()}
    results = pd.DataFrame([(rnA1[r], rnA2[c], s) for r, c, s in aligned],
                           columns=[a, b, "score"])

    tsv = os.path.join(RESULTS, f"{a}_{b}_alignment_scored.tsv")
    results.to_csv(tsv, sep="\t", index=False)
    smaller_n = min(len(nA1), len(nA2))
    print(f"  wrote {tsv} ({len(results)} pairs, "
          f"{100 * len(results) / smaller_n:.1f}% of the smaller network)")

    scores = results["score"].to_numpy()
    summary = {
        "species_a": a, "species_b": b,
        "n_pairs": len(scores), "n_a_nodes": len(nA1), "n_b_nodes": len(nA2),
        "coverage_of_smaller": len(scores) / smaller_n,
        "alpha": ALPHA, "niter": NITER,
        "score_min": float(scores.min()), "score_max": float(scores.max()),
        "score_median": float(np.median(scores)),
        "score_p10": float(np.percentile(scores, 10)),
        "score_p20": float(np.percentile(scores, 20)),
    }
    js = os.path.join(RESULTS, f"{a}_{b}_alignment_scores.json")
    with open(js, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  wrote {js}")


if __name__ == "__main__":
    main()
