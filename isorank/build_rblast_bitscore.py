#!/usr/bin/env python3
"""
Same merge as build_rblast.py (union of both BLAST directions, higher
value kept per pair, then min-max rescaled to [0, 1]) but using raw
bitscore directly as the similarity signal instead of -log10(evalue).
Written to compare against the evalue-based score used elsewhere in this
project -- bitscore is monotonic in evalue for a fixed alignment length but
isn't identically distributed, so this checks whether the choice of raw
signal changes the alignment's behavior.

Usage:
    python3 build_rblast_bitscore.py <a> <b> <alias_a> <alias_b>
Expects blast_out/<a>_vs_<b>.tsv and blast_out/<b>_vs_<a>.tsv (tabular
outfmt 6: qseqid sseqid bitscore evalue pident).
Writes net/<alias_a>-<alias_b>.tsv with header alias_a/alias_b/score.
"""

import sys
import pandas as pd

a, b, alias_a, alias_b = sys.argv[1:5]
cols = ["qseqid", "sseqid", "bitscore", "evalue", "pident"]

fwd = pd.read_csv(f"blast_out/{a}_vs_{b}.tsv", sep="\t", names=cols)
fwd = fwd.rename(columns={"qseqid": a, "sseqid": b})[[a, b, "bitscore"]]

rev = pd.read_csv(f"blast_out/{b}_vs_{a}.tsv", sep="\t", names=cols)
rev = rev.rename(columns={"qseqid": b, "sseqid": a})[[a, b, "bitscore"]]

both = pd.concat([fwd, rev], ignore_index=True)

merged = both.groupby([a, b], as_index=False)["bitscore"].max()
lo, hi = merged["bitscore"].min(), merged["bitscore"].max()
merged["score"] = (merged["bitscore"] - lo) / (hi - lo)
merged = merged[[a, b]].join(merged["score"])
merged = merged.rename(columns={a: alias_a, b: alias_b})

merged.to_csv(f"net/{alias_a}-{alias_b}.tsv", sep="\t", index=False)
print(f"forward hits: {len(fwd)}, reverse hits: {len(rev)}, "
      f"merged unique pairs: {len(merged)}, raw bitscore range [{lo:.2f}, {hi:.2f}]")
