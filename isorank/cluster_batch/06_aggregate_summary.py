#!/usr/bin/env python3
"""
Collect every <a>_lcc_bs_<b>_lcc_bs_alignment_scores.json produced by
05_run_all_pairs_alignment.sbatch (up to 45, one per species pair) into a
single summary table -- run this after the batch alignment job (it's also
called automatically at the end of that sbatch script).

Usage:
    python3 06_aggregate_summary.py
Expects species.tsv and ../<a>_lcc_bs_<b>_lcc_bs_alignment_scores.json for each pair.
Writes ../all_pairs_summary.tsv and ../all_pairs_summary.json
"""

import csv
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ISORANK_DIR = os.path.dirname(HERE)

rows = []
for path in sorted(glob.glob(os.path.join(ISORANK_DIR, "*_lcc_bs_*_alignment_scores.json"))):
    with open(path) as f:
        d = json.load(f)
    rows.append(d)

missing = []
shorts = []
with open(os.path.join(HERE, "species.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    shorts = [r["short"] for r in reader]

expected = {tuple(sorted((shorts[i], shorts[j])))
            for i in range(len(shorts)) for j in range(i + 1, len(shorts))}
found = {tuple(sorted((r["species_a"].removesuffix("_lcc_bs"), r["species_b"].removesuffix("_lcc_bs"))))
         for r in rows}
missing = sorted(expected - found)

print(f"found {len(rows)}/{len(expected)} pair results")
if missing:
    print(f"missing {len(missing)} pairs: {missing}")

out_tsv = os.path.join(ISORANK_DIR, "all_pairs_summary.tsv")
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species_a", "species_b", "n_pairs", "n_a_nodes", "n_b_nodes",
                "score_min", "score_max", "score_median", "score_p10", "score_p20"])
    for r in rows:
        w.writerow([r["species_a"].removesuffix("_lcc_bs"), r["species_b"].removesuffix("_lcc_bs"),
                    r["n_pairs"], r["n_a_nodes"], r["n_b_nodes"],
                    r["score_min"], r["score_max"], r["score_median"],
                    r.get("score_p10"), r.get("score_p20")])
print(f"wrote {out_tsv}")

out_json = os.path.join(ISORANK_DIR, "all_pairs_summary.json")
with open(out_json, "w") as f:
    json.dump({"n_found": len(rows), "n_expected": len(expected),
               "missing_pairs": missing, "results": rows}, f, indent=2)
print(f"wrote {out_json}")
