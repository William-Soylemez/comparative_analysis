#!/usr/bin/env python3
"""Merge reciprocal BLASTP hits (species_a<->species_b) into netalign's expected
3-column rblast format: header row `species_a\tspecies_b\tscore`, one row per
(species_a_protein, species_b_protein) pair.

Score follows the original IsoRank paper's convention (not raw bitscore):
    raw   = -log10(evalue), evalue floored at 1e-300 to avoid -log10(0) on
            BLAST's underflowed-to-zero e-values for extremely strong hits
    score = raw min-max rescaled to [0, 1] across all merged pairs
If a pair was hit in both BLAST directions, the higher (more significant,
i.e. lower e-value) raw score is kept before rescaling.

Usage:
    python3 build_rblast.py <species_a> <species_b>
Expects blast_out/<species_a>_vs_<species_b>.tsv and
        blast_out/<species_b>_vs_<species_a>.tsv (tabular outfmt 6:
        qseqid sseqid bitscore evalue pident)
Writes net/<species_a>-<species_b>.tsv
"""

import sys
import numpy as np
import pandas as pd

EVALUE_FLOOR = 1e-300

a, b = sys.argv[1], sys.argv[2]
cols = ["qseqid", "sseqid", "bitscore", "evalue", "pident"]

fwd = pd.read_csv(f"blast_out/{a}_vs_{b}.tsv", sep="\t", names=cols)
fwd = fwd.rename(columns={"qseqid": a, "sseqid": b})[[a, b, "evalue"]]

rev = pd.read_csv(f"blast_out/{b}_vs_{a}.tsv", sep="\t", names=cols)
rev = rev.rename(columns={"qseqid": b, "sseqid": a})[[a, b, "evalue"]]

both = pd.concat([fwd, rev], ignore_index=True)
both["raw"] = -np.log10(both["evalue"].clip(lower=EVALUE_FLOOR))

merged = both.groupby([a, b], as_index=False)["raw"].max()
lo, hi = merged["raw"].min(), merged["raw"].max()
merged["score"] = (merged["raw"] - lo) / (hi - lo)
merged = merged[[a, b, "score"]]

merged.to_csv(f"net/{a}-{b}.tsv", sep="\t", index=False)
print(f"forward hits: {len(fwd)}, reverse hits: {len(rev)}, "
      f"merged unique pairs: {len(merged)}, raw -log10(evalue) range [{lo:.2f}, {hi:.2f}]")
