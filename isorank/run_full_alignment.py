#!/usr/bin/env python3
"""
Full IsoRank alignment between two species, matching 100% of the smaller
species' nodes and recording the IsoRank similarity score for every matched
pair.

Calls netalign's own compute_adjacency / compute_pairs / compute_isorank
directly (rather than the `netalign isorank` CLI) because the CLI's greedy
assignment only returns index pairs, not the score at each match — we need
the score per pair to later choose a confidence cutoff.

Usage:
    python3 run_full_alignment.py <species_a> <species_b>
Expects net/<species_a>.tsv, net/<species_b>.tsv, net/<species_a>-<species_b>.tsv
Outputs:
    <species_a>_<species_b>_alignment_scored.tsv (one row per node of the
        smaller network, i.e. 100% coverage of the smaller species)
    <species_a>_<species_b>_alignment_scores.json (summary stats)
"""

import sys
import json
import numpy as np
import pandas as pd

from netalign.approx_isorank.io_utils import compute_adjacency, compute_pairs
from netalign.approx_isorank.isorank_compute import compute_isorank

ALPHA = 0.7
NITER = 20

a, b = sys.argv[1], sys.argv[2]

print(f"Loading networks and reciprocal-BLAST pairs ({a} <-> {b}) ...")
df1 = pd.read_csv(f"net/{a}.tsv", sep="\t", header=None)
df2 = pd.read_csv(f"net/{b}.tsv", sep="\t", header=None)
dpairs = pd.read_csv(f"net/{a}-{b}.tsv", sep="\t")

Af1, nA1 = compute_adjacency(df1)
Af2, nA2 = compute_adjacency(df2)
print(f"  {a}: {len(nA1)} nodes, {b}: {len(nA2)} nodes")

E = compute_pairs(dpairs, nA1, nA2, a, b)

print("Computing Isorank similarity matrix ...")
R = compute_isorank(Af1, Af2, E, alpha=ALPHA, maxiter=NITER)[-1]
print(f"  R shape: {R.shape}")

print("Greedy one-to-one assignment, recording score per pair, "
      "matching 100% of the smaller network ...")
n_align = min(R.shape)  # = size of the smaller species -> full coverage
Rwork = R.copy()
aligned = []  # (row_idx, col_idx, score)
for _ in range(n_align):
    maxcols = np.argmax(Rwork, axis=1)
    row_best_vals = Rwork[np.arange(Rwork.shape[0]), maxcols]
    maxid = int(np.argmax(row_best_vals))
    maxcol = int(maxcols[maxid])
    score = float(Rwork[maxid, maxcol])
    aligned.append((maxid, maxcol, score))
    Rwork[:, maxcol] = -1
    Rwork[maxid, :] = -1

rnA1 = {v: k for k, v in nA1.items()}
rnA2 = {v: k for k, v in nA2.items()}

results = pd.DataFrame(
    [(rnA1[r], rnA2[c], s) for r, c, s in aligned],
    columns=[a, b, "score"],
)
outfile = f"{a}_{b}_alignment_scored.tsv"
results.to_csv(outfile, sep="\t", index=False)
smaller_n = min(len(nA1), len(nA2))
print(f"  wrote {outfile} ({len(results)} pairs, "
      f"{100*len(results)/smaller_n:.1f}% of the smaller species' nodes)")

scores = results["score"].to_numpy()
summary = {
    "species_a": a,
    "species_b": b,
    "n_pairs": len(scores),
    "n_a_nodes": len(nA1),
    "n_b_nodes": len(nA2),
    "coverage_of_smaller_species": len(scores) / smaller_n,
    "alpha": ALPHA,
    "niter": NITER,
    "score_min": float(scores.min()),
    "score_max": float(scores.max()),
    "score_median": float(np.median(scores)),
    "score_p10": float(np.percentile(scores, 10)),
    "score_p20": float(np.percentile(scores, 20)),
}
jsonfile = f"{a}_{b}_alignment_scores.json"
with open(jsonfile, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  wrote {jsonfile}")
print(json.dumps(summary, indent=2))
