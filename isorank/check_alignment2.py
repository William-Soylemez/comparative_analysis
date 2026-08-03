#!/usr/bin/env python3
"""Compare IsoRank's top 2000 pairs against a pure-BLAST top-2000 baseline
(no network topology at all — just best reciprocal-BLAST hit per protein,
ranked by score), on GO-term sharing. Reports both:
  - RAW: any shared GO term from GO_map.csv as-is (includes generic/redundant
    parent terms, e.g. "metabolic process").
  - LEAF: same check but each protein's annotation is first pruned to only
    its most-specific ("leaf") terms per species, using the same is_a/part_of
    ancestor-pruning as analysis_c3_leaf_terms.py — removes the redundant
    generic parents so a shared hit means something more specific.
"""

import sys
import csv
import random

sys.path.insert(0, "..")
from goslim_util import _parse_obo

random.seed(0)

OBO = "../go_data/go-basic.obo"


class Ontology:
    def __init__(self, obo_path):
        self.name, self.namespace, self._parents, self.obsolete, self._alt = _parse_obo(obo_path)
        self._anc_cache = {}

    def resolve(self, term):
        return self._alt.get(term, term)

    def ancestors(self, term):
        term = self.resolve(term)
        if term in self._anc_cache:
            return self._anc_cache[term]
        seen, stack = set(), [term]
        while stack:
            t = stack.pop()
            for p in self._parents.get(t, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        self._anc_cache[term] = seen
        return seen


def leaf_terms(all_terms, onto):
    redundant = set()
    for t in all_terms:
        redundant |= (onto.ancestors(t) & all_terms)
    return all_terms - redundant


def load_go(path):
    p2terms = {}
    all_terms = set()
    with open(path) as f:
        for row in csv.DictReader(f):
            terms = {g for g in (row.get("GO_list") or "").split(";") if g}
            p2terms[row["prot_id"]] = terms
            all_terms |= terms
    return p2terms, all_terms


print("Loading go-basic.obo ...")
onto = Ontology(OBO)

go_n_raw, all_n = load_go("../input/GCF_000182925.2/GCF_000182925.2_GO_map.csv")
go_c_raw, all_c = load_go("../input/GCF_000182965.3/GCF_000182965.3_GO_map.csv")

leaves_n = leaf_terms(all_n, onto)
leaves_c = leaf_terms(all_c, onto)
print(f"  ncrassa: {len(all_n)} annotated terms -> {len(leaves_n)} leaf terms")
print(f"  calbicans: {len(all_c)} annotated terms -> {len(leaves_c)} leaf terms")

go_n_leaf = {p: (terms & leaves_n) for p, terms in go_n_raw.items()}
go_c_leaf = {p: (terms & leaves_c) for p, terms in go_c_raw.items()}


def share(p2go_a, p2go_b, n, c):
    return len(p2go_a.get(n, set()) & p2go_b.get(c, set())) > 0


def pct_share(pairs, p2go_a, p2go_b):
    hits = sum(1 for n, c in pairs if share(p2go_a, p2go_b, n, c))
    return hits, len(pairs), 100 * hits / len(pairs)


def random_baseline(pairs, p2go_a, p2go_b, n_shuffles=20):
    c_pool = [c for _, c in pairs]
    vals = []
    for _ in range(n_shuffles):
        shuffled = c_pool[:]
        random.shuffle(shuffled)
        hits = sum(1 for (n, _), c in zip(pairs, shuffled) if share(p2go_a, p2go_b, n, c))
        vals.append(hits)
    avg = sum(vals) / len(vals)
    return avg, 100 * avg / len(pairs)


# --- IsoRank pairs ---
isorank_pairs = []
with open("results/ncrassa_calbicans_alignment.tsv") as f:
    r = csv.reader(f, delimiter="\t")
    next(r)
    for row in r:
        isorank_pairs.append((row[0], row[1]))

# --- Pure-BLAST top-2000 baseline (no topology): best hit per ncrassa protein, ranked by score ---
best = {}
with open("net/ncrassa-calbicans.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        n, c, s = row["ncrassa"], row["calbicans"], float(row["score"])
        if n not in best or s > best[n][1]:
            best[n] = (c, s)

ranked = sorted(best.items(), key=lambda kv: kv[1][1], reverse=True)
blast_pairs = [(n, c) for n, (c, s) in ranked[:len(isorank_pairs)]]

print(f"\nN pairs compared: {len(isorank_pairs)} (IsoRank) vs {len(blast_pairs)} (pure-BLAST top-N)\n")

for label, pairs in [("IsoRank (topology + BLAST)", isorank_pairs),
                     ("Pure BLAST top-N (no topology)", blast_pairs)]:
    print(f"=== {label} ===")
    h, n, p = pct_share(pairs, go_n_raw, go_c_raw)
    print(f"  RAW GO sharing:  {h}/{n} ({p:.1f}%)")
    h, n, p = pct_share(pairs, go_n_leaf, go_c_leaf)
    print(f"  LEAF GO sharing: {h}/{n} ({p:.1f}%)")
    avg, p = random_baseline(pairs, go_n_raw, go_c_raw)
    print(f"  random baseline (raw, 20 shuffles): {avg:.1f}/{n} ({p:.1f}%)")
    avg, p = random_baseline(pairs, go_n_leaf, go_c_leaf)
    print(f"  random baseline (leaf, 20 shuffles): {avg:.1f}/{n} ({p:.1f}%)")
    print()
