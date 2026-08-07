#!/usr/bin/env python3
"""
Step 2/3: turn reciprocal BLAST into the two inputs IsoRank needs for a pair.

For each species:
  * read its PHILHARMONIC network (`<acc>_network.positive.tsv`), keep only the
    two protein-id columns (the 3rd is a D-SCRIPT score IsoRank doesn't use),
    restrict to the largest connected component, and write net/<short>.tsv.

For the pair:
  * merge the two BLAST directions, keep the max bitscore per protein pair, and
    min-max normalize it to [0, 1] -- this normalized bitscore is the sequence-
    similarity signal IsoRank blends with network topology. Write net/<a>-<b>.tsv.

Usage: python3 2_similarity.py <a> <b>
  e.g. python3 2_similarity.py scer calbicans
Expects: blast_out/<a>_vs_<b>.tsv, blast_out/<b>_vs_<a>.tsv (from 1_blast.sh),
         species.txt (short -> accession),
         ../input/<acc>/<acc>_network.positive.tsv
Writes:  net/<a>.tsv, net/<b>.tsv (LCC edge lists), net/<a>-<b>.tsv (similarity)
"""

import csv
import os
import sys

import pandas as pd
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
NET = os.path.join(HERE, "net")
INPUT = os.path.join(os.path.dirname(HERE), "input")


def accession(short):
    with open(os.path.join(HERE, "species.txt")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["short"] == short:
                return row["accession"]
    sys.exit(f"{short}: not found in species.txt")


def write_lcc_edges(short):
    """Extract the 2-column edge list from the species' network, restrict to the
    largest connected component, and write net/<short>.tsv."""
    acc = accession(short)
    net_path = os.path.join(INPUT, acc, f"{acc}_network.positive.tsv")
    if not os.path.exists(net_path):
        sys.exit(f"missing network file: {net_path}")

    g = nx.Graph()
    with open(net_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != p[1]:
                g.add_edge(p[0], p[1])

    lcc = max(nx.connected_components(g), key=len)
    sub = g.subgraph(lcc)
    out = os.path.join(NET, f"{short}.tsv")
    with open(out, "w") as f:
        for u, v in sub.edges():
            f.write(f"{u}\t{v}\n")
    print(f"{short}: {g.number_of_nodes()} nodes / {g.number_of_edges()} edges -> "
          f"LCC {sub.number_of_nodes()} nodes / {sub.number_of_edges()} edges -> {out}")


def build_similarity(a, b):
    """Merge both BLAST directions -> max bitscore per pair -> normalized [0,1]."""
    cols = ["qseqid", "sseqid", "bitscore", "evalue", "pident"]

    def hits(q, s, qname, sname):
        df = pd.read_csv(os.path.join(HERE, "blast_out", f"{q}_vs_{s}.tsv"), sep="\t", names=cols)
        return df.rename(columns={"qseqid": qname, "sseqid": sname})[[qname, sname, "bitscore"]]

    both = pd.concat([hits(a, b, a, b), hits(b, a, b, a)], ignore_index=True)
    merged = both.groupby([a, b], as_index=False)["bitscore"].max()
    lo, hi = merged["bitscore"].min(), merged["bitscore"].max()
    merged["score"] = (merged["bitscore"] - lo) / (hi - lo)

    out = os.path.join(NET, f"{a}-{b}.tsv")
    merged[[a, b, "score"]].to_csv(out, sep="\t", index=False)
    print(f"similarity: {len(merged)} unique pairs, raw bitscore [{lo:.1f}, {hi:.1f}] -> {out}")


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: python3 2_similarity.py <a> <b>")
    a, b = sys.argv[1], sys.argv[2]
    os.makedirs(NET, exist_ok=True)
    write_lcc_edges(a)
    write_lcc_edges(b)
    build_similarity(a, b)


if __name__ == "__main__":
    main()
