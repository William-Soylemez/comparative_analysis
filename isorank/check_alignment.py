#!/usr/bin/env python3
"""Quick sanity-check of the IsoRank alignment: (1) how many pairs match the
best BLAST hit alone (topology's added value), (2) GO-term overlap enrichment
vs a random baseline."""

import csv
import random
from collections import defaultdict

random.seed(0)

align = []
with open("ncrassa_calbicans_alignment.tsv") as f:
    r = csv.reader(f, delimiter="\t")
    next(r)
    for row in r:
        align.append((row[0], row[1]))

# best BLAST hit per ncrassa protein (max bitscore), from the merged rblast file
best_hit = {}
with open("net/ncrassa-calbicans.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    best = {}
    for row in r:
        n, c, s = row["ncrassa"], row["calbicans"], float(row["score"])
        if n not in best or s > best[n][1]:
            best[n] = (c, s)
    best_hit = {n: c for n, (c, s) in best.items()}

matches_best_blast = sum(1 for n, c in align if best_hit.get(n) == c)
print(f"IsoRank pairs: {len(align)}")
print(f"Pairs that ARE also the top BLAST hit: {matches_best_blast} "
      f"({100*matches_best_blast/len(align):.1f}%)")
print(f"Pairs where topology picked a DIFFERENT partner than best BLAST hit: "
      f"{len(align)-matches_best_blast} ({100*(len(align)-matches_best_blast)/len(align):.1f}%)")

# GO overlap check
def load_go(path):
    p2go = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            p2go[row["prot_id"]] = {g for g in (row.get("GO_list") or "").split(";") if g}
    return p2go

go_n = load_go("../GCF_000182925.2/GCF_000182925.2_GO_map.csv")
go_c = load_go("../GCF_000182965.3/GCF_000182965.3_GO_map.csv")

def share_go(n, c):
    return len(go_n.get(n, set()) & go_c.get(c, set())) > 0

aligned_share = sum(1 for n, c in align if share_go(n, c))
print(f"\nAligned pairs sharing >=1 GO term: {aligned_share}/{len(align)} "
      f"({100*aligned_share/len(align):.1f}%)")

# random baseline: shuffle the calbicans side
c_pool = [c for _, c in align]
random_shares = []
for _ in range(20):
    shuffled = c_pool[:]
    random.shuffle(shuffled)
    s = sum(1 for (n, _), c in zip(align, shuffled) if share_go(n, c))
    random_shares.append(s)
avg_random = sum(random_shares) / len(random_shares)
print(f"Random-pair baseline sharing >=1 GO term: {avg_random:.1f}/{len(align)} "
      f"({100*avg_random/len(align):.1f}%, avg of 20 shuffles)")
