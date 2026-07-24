#!/usr/bin/env python3
"""
Sanity check: run the netalign package's own logic (same functions the
approx_isorank_run_matching.ipynb notebook calls) on the notebook's own real
example data (bakers.s.tsv / rat.s.tsv / rat-bakers.tsv from the repo's
data/intact/), to see what the final score distribution looks like on data
the paper authors themselves published — for direct comparison against our
own species pairs, which topped out around 2.5e-4.

Mirrors run_full_alignment.py: full one-to-one coverage of the smaller
network (bakers, 6478 nodes), alpha=0.7, niter=20 (the "close to true
IsoRank" setting per their docs, not the niter=5 demo in the notebook).
"""

import json
import numpy as np
import pandas as pd

from netalign.approx_isorank.io_utils import compute_adjacency, compute_pairs
from netalign.approx_isorank.isorank_compute import compute_isorank

DATA = "/tmp/netalign_repo_check/data/intact"
ALPHA = 0.7
NITER = 20

print("Loading bakers/rat networks + rat-bakers reciprocal-BLAST pairs ...")
df1 = pd.read_csv(f"{DATA}/bakers.s.tsv", sep="\t", header=None)
df2 = pd.read_csv(f"{DATA}/rat.s.tsv", sep="\t", header=None)
dpairs = pd.read_csv(f"{DATA}/rat-bakers.tsv", sep="\t")

org1, org2 = "bakers", "rat"

Af1, nA1 = compute_adjacency(df1)
Af2, nA2 = compute_adjacency(df2)
print(f"  bakers: {len(nA1)} nodes, rat: {len(nA2)} nodes")

E = compute_pairs(dpairs, nA1, nA2, org1, org2)

print("Computing Isorank similarity matrix ...")
R = compute_isorank(Af1, Af2, E, alpha=ALPHA, maxiter=NITER)[-1]
print(f"  R shape: {R.shape}")

print("Greedy one-to-one assignment, recording score per pair, "
      "matching 100% of the smaller network ...")
n_align = min(R.shape)
Rwork = R.copy()
aligned = []
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
    columns=[org1, org2, "score"],
)
results.to_csv("bakers_rat_alignment_scored.tsv", sep="\t", index=False)
print(f"  wrote bakers_rat_alignment_scored.tsv ({len(results)} pairs)")

scores = results["score"].to_numpy()
summary = {
    "species_a": org1, "species_b": org2,
    "n_pairs": len(scores), "n_a_nodes": len(nA1), "n_b_nodes": len(nA2),
    "alpha": ALPHA, "niter": NITER,
    "score_min": float(scores.min()), "score_max": float(scores.max()),
    "score_median": float(np.median(scores)),
    "score_p90": float(np.percentile(scores, 90)),
    "score_p99": float(np.percentile(scores, 99)),
}
with open("bakers_rat_alignment_scores.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
