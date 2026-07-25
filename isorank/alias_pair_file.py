#!/usr/bin/env python3
"""
Copy an existing net/<a>-<b>.tsv reciprocal-BLAST pair file to
net/<alias_a>-<alias_b>.tsv with its header columns renamed, so
run_full_alignment.py (which reads net/<a>.tsv, net/<b>.tsv, net/<a>-<b>.tsv
purely by name) can be pointed at an aliased network (e.g. "<sp>_lcc") while
reusing the same BLAST-derived pairs unchanged. compute_pairs() already
drops any pair referencing a node absent from the aliased network's
nodemap, so no row filtering is needed here -- just the rename.

Usage:
    python3 alias_pair_file.py <a> <b> <alias_a> <alias_b>
Expects net/<a>-<b>.tsv. Writes net/<alias_a>-<alias_b>.tsv.
"""

import sys
import pandas as pd

a, b, alias_a, alias_b = sys.argv[1:5]

df = pd.read_csv(f"net/{a}-{b}.tsv", sep="\t")
df = df.rename(columns={a: alias_a, b: alias_b})
df.to_csv(f"net/{alias_a}-{alias_b}.tsv", sep="\t", index=False)
print(f"wrote net/{alias_a}-{alias_b}.tsv ({len(df)} pairs)")
