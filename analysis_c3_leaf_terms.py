#!/usr/bin/env python3
"""
Analysis C3 — cluster spread of non-redundant ("leaf") GO terms.

Supervisor-suggested alternative to Analysis C's concentration score:
  1. For each species, take the set of GO terms actually annotated to at least
     one protein. Drop any term that is an ANCESTOR (via is_a + part_of,
     transitively) of another annotated term — i.e. keep only the most
     specific ("leaf") term along each annotated lineage. This removes the
     redundancy where a generic parent term (e.g. "metabolic process") is
     annotated purely because a more specific child term (e.g. "chitin
     catabolic process") is annotated to the same proteins.
  2. For each surviving leaf term, count how many DISTINCT clusters it
     appears in (among clustered proteins only).
  3. Plot the distribution of that per-term cluster count, per species.

This is a simpler, more direct read than Analysis C's frac_biggest: instead of
"how concentrated is this term's protein set", it asks "how many modules does
this term's protein set touch, after removing terms that are just redundant
restatements of a more specific child term."

Usage:
    python3 analysis_c3_leaf_terms.py [SPECIES_DIR ...]

Outputs go_leaf_spread.json, go_leaf_spread.png.
"""

import sys
import os
import re
import csv
import json
import glob
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from goslim_util import _parse_obo

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]
OBO = "go_data/go-basic.obo"
MIN_PROTEINS = 8  # same floor as Analysis C — drop tiny terms where n_clusters is trivially small


def discover_species(args):
    dirs = args if args else sorted(
        d for d in os.listdir(".") if os.path.isdir(d) and glob.glob(f"{d}/*_clusters.json")
    )
    out = []
    for d in dirs:
        d = d.rstrip("/")
        hits = glob.glob(f"{d}/*_clusters.json")
        if hits:
            out.append((os.path.basename(hits[0]).replace("_clusters.json", ""), d))
    return out


def scrape_go_names(species):
    names = {}
    pat = re.compile(r"(GO:\d+)\s*-\s*<([^>]+)>")
    for _, d in species:
        for hr in glob.glob(f"{d}/*_human_readable.txt"):
            with open(hr) as f:
                for m in pat.finditer(f.read()):
                    names.setdefault(m.group(1), m.group(2))
    return names


def protein_to_cluster(clusters_path):
    p2c = {}
    for cid, c in json.load(open(clusters_path)).items():
        for m in c["members"]:
            p2c[m] = cid
    return p2c


class Ontology:
    """Thin wrapper around goslim_util's OBO parser for ancestor lookups
    (no slim subset involved — this is the full ontology)."""

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


def load_go_map(path):
    """Return protein -> set of resolved GO terms, and the union of all terms seen."""
    p2terms = {}
    all_terms = set()
    with open(path) as f:
        for row in csv.DictReader(f):
            terms = {g for g in (row.get("GO_list") or "").split(";") if g}
            p2terms[row["prot_id"]] = terms
            all_terms |= terms
    return p2terms, all_terms


def leaf_terms(all_terms, onto):
    """Drop any term that is an ancestor of another term in the same set."""
    redundant = set()
    for t in all_terms:
        redundant |= (onto.ancestors(t) & all_terms)
    return all_terms - redundant, redundant


def score_species(name, d, onto):
    p2c = protein_to_cluster(f"{d}/{name}_clusters.json")
    p2terms, all_terms = load_go_map(f"{d}/{name}_GO_map.csv")

    leaves, dropped = leaf_terms(all_terms, onto)

    term_clusters = defaultdict(set)
    term_nproteins = Counter()
    for pid, terms in p2terms.items():
        cid = p2c.get(pid)
        if cid is None:
            continue
        for g in terms:
            if g not in leaves:
                continue
            term_clusters[g].add(cid)
            term_nproteins[g] += 1

    scores = {}
    for g, cids in term_clusters.items():
        n = term_nproteins[g]
        if n < MIN_PROTEINS:
            continue
        scores[g] = {"n_proteins": n, "n_clusters": len(cids)}
    return scores, len(all_terms), len(leaves), len(dropped)


def plot(all_scores, species_names, outfile):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    for i, sp in enumerate(species_names):
        color = PALETTE[i % len(PALETTE)]
        counts = np.array([v["n_clusters"] for v in all_scores[sp].values()])

        ax1.hist(counts, bins=range(1, counts.max() + 2), histtype="step",
                 linewidth=2, color=color, label=sp)

        vals = np.sort(counts)
        ccdf = 1.0 - np.arange(len(vals)) / len(vals)
        ax2.plot(vals, ccdf, linewidth=2, color=color, label=sp)

    ax1.set_yscale("log")
    ax1.set_xlabel("Number of distinct clusters a leaf GO term appears in")
    ax1.set_ylabel("Number of GO terms (log)")
    ax1.set_title("Cluster-spread histogram")
    ax1.grid(True, which="both", color="#e6e6e3", linewidth=0.6)
    ax1.set_axisbelow(True)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Clusters k (log)")
    ax2.set_ylabel("Fraction of leaf terms spanning ≥ k clusters")
    ax2.set_title("Cluster-spread (CCDF, size-fair)")
    ax2.grid(True, which="both", color="#e6e6e3", linewidth=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, title="Species")

    fig.suptitle(f"Analysis C3 — how many clusters does a non-redundant GO term touch?\n"
                 f"(parent terms of another annotated term dropped; ≥{MIN_PROTEINS} proteins)",
                 fontsize=12, y=1.0)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No clusters.json files found.", file=sys.stderr)
        sys.exit(1)
    print("Loading go-basic.obo ...")
    onto = Ontology(OBO)
    go_names = scrape_go_names(species)

    species_names = [s[0] for s in species]
    all_scores = {}
    for name, d in species:
        sc, n_all, n_leaf, n_dropped = score_species(name, d, onto)
        all_scores[name] = sc
        counts = [v["n_clusters"] for v in sc.values()]
        print(f"  {name}: {n_all} annotated terms -> {n_leaf} leaf terms "
              f"({n_dropped} redundant parents dropped); {len(sc)} scored (>={MIN_PROTEINS} proteins); "
              f"median clusters/term {np.median(counts):.1f}, max {max(counts)}")

    plot(all_scores, species_names, "go_leaf_spread.png")

    out = {
        "params": {"min_proteins": MIN_PROTEINS, "obo": OBO},
        "summary": {
            sp: {
                "n_terms_scored": len(all_scores[sp]),
                "median_n_clusters": float(np.median([v["n_clusters"] for v in all_scores[sp].values()])),
                "max_n_clusters": max(v["n_clusters"] for v in all_scores[sp].values()),
            }
            for sp in species_names
        },
        "per_species_terms": {
            sp: {g: {**v, "name": go_names.get(g, "")} for g, v in all_scores[sp].items()}
            for sp in species_names
        },
    }
    with open("go_leaf_spread.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote go_leaf_spread.json")


if __name__ == "__main__":
    main()
