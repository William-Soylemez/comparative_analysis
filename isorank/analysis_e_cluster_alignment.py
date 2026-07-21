#!/usr/bin/env python3
"""
Analysis E — cluster alignment: is a species-A module's internal wiring
conserved in species B?

For each cluster C_A in species A:
  1. Map its members through the IsoRank protein alignment -> corresponding
     nodes C_B in species B (only members with an aligned partner count).
  2. Induce G_ca = the subgraph of species A's network restricted to the
     mapped subset of C_A (induced subgraph = same nodes, only edges where
     BOTH endpoints are in the set).
  3. Induce G_cb = the subgraph of species B's network restricted to C_B.
  4. A species-A edge (u, v) is "conserved" if (map(u), map(v)) is also an
     edge in G_cb.
  5. Report three standard network-alignment edge-conservation scores:
       EC  (edge correctness)          = conserved / |E(G_ca)|
       ICS (induced conserved structure) = conserved / |E(G_cb)|
       S3  (symmetric substructure score) = conserved / (|E(G_ca)| + |E(G_cb)| - conserved)
     S3 is the most balanced (doesn't privilege either species' density) and
     is the headline number; EC/ICS are reported to show if conservation is
     lopsided.

Clusters need >= MIN_MAPPED proteins with an alignment partner, AND both
induced subgraphs need >= 1 internal edge, to be scored (small/edgeless
clusters give meaningless 0/1 or 1/1 ratios).

Usage:
    python3 analysis_e_cluster_alignment.py
        [--net1_dir DIR] [--net2_dir DIR] [--alignment FILE] [--min_mapped N]

Defaults assume this script lives in isorank/ next to net/ and the alignment
output, with species dirs one level up.

Outputs cluster_alignment.json + cluster_alignment.png.
"""

import argparse
import csv
import json
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]


def load_alignment(path):
    m = {}
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        for row in r:
            m[row[0]] = row[1]
    return m


def load_clusters(path):
    return {cid: c["members"] for cid, c in json.load(open(path)).items()}


def load_edges(path):
    """Two-or-more column tab-delimited network file -> dict node -> set(neighbors)."""
    adj = defaultdict(set)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            if a == b:
                continue
            adj[a].add(b)
            adj[b].add(a)
    return adj


def induced_edges(nodes, adj):
    """Edges (as frozensets) among `nodes`, using full-network adjacency `adj`."""
    nodes = set(nodes)
    edges = set()
    for u in nodes:
        for v in adj.get(u, ()):
            if v in nodes and u < v:
                edges.add((u, v))
            elif v in nodes and v < u:
                edges.add((v, u))
    return edges


def score_cluster(members_a, align, adj_a, adj_b, min_mapped):
    mapped = [(u, align[u]) for u in members_a if u in align]
    if len(mapped) < min_mapped:
        return None
    c_a = {u for u, v in mapped}
    c_b = {v for u, v in mapped}
    u2v = dict(mapped)

    e_ca = induced_edges(c_a, adj_a)
    e_cb = induced_edges(c_b, adj_b)
    if not e_ca or not e_cb:
        return None

    conserved = 0
    for u, v in e_ca:
        mu, mv = u2v[u], u2v[v]
        if mv in adj_b.get(mu, ()):
            conserved += 1

    ec = conserved / len(e_ca)
    ics = conserved / len(e_cb)
    denom = len(e_ca) + len(e_cb) - conserved
    s3 = conserved / denom if denom > 0 else 0.0

    return {
        "n_members": len(members_a),
        "n_mapped": len(mapped),
        "n_edges_ca": len(e_ca),
        "n_edges_cb": len(e_cb),
        "conserved_edges": conserved,
        "EC": ec,
        "ICS": ics,
        "S3": s3,
    }


def plot(scores, outfile, label_a, label_b, min_mapped):
    metrics = ["EC", "ICS", "S3"]
    data = [[s[m] for s in scores] for m in metrics]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(PALETTE[i % len(PALETTE)])
        body.set_edgecolor(PALETTE[i % len(PALETTE)])
        body.set_alpha(0.55)
    parts["cmedians"].set_color("#333333")
    for i, vals in enumerate(data):
        ax.text(i + 1, np.median(vals), f"  med {np.median(vals):.2f}",
                va="center", fontsize=9, color="#333333")
        ax.scatter(np.random.normal(i + 1, 0.04, len(vals)), vals,
                   s=10, color=PALETTE[i % len(PALETTE)], alpha=0.4, zorder=3)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["EC\n(vs A's edges)", "ICS\n(vs B's edges)", "S3\n(symmetric)"])
    ax.set_ylabel("Edge conservation score")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"Analysis E — cluster wiring conservation, {label_a} -> {label_b}\n"
                 f"({len(scores)} clusters scored, >={min_mapped} mapped proteins & >=1 internal edge each side)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net1_dir", default="../GCF_000182925.2")
    ap.add_argument("--net2_dir", default="../GCF_000182965.3")
    ap.add_argument("--net1_name", default="GCF_000182925.2")
    ap.add_argument("--net2_name", default="GCF_000182965.3")
    ap.add_argument("--alignment", default="ncrassa_calbicans_alignment_full.tsv")
    ap.add_argument("--min_mapped", type=int, default=5)
    args = ap.parse_args()

    print(f"Loading alignment ({args.alignment}) ...")
    align = load_alignment(args.alignment)
    print(f"  {len(align)} aligned protein pairs")

    clusters_a = load_clusters(f"{args.net1_dir}/{args.net1_name}_clusters.json")
    adj_a = load_edges(f"{args.net1_dir}/{args.net1_name}_network.positive.tsv")
    adj_b = load_edges(f"{args.net2_dir}/{args.net2_name}_network.positive.tsv")
    print(f"  {len(clusters_a)} clusters in {args.net1_name}")

    scores = []
    n_skipped_coverage = 0
    n_skipped_no_edges = 0
    for cid, members in clusters_a.items():
        s = score_cluster(members, align, adj_a, adj_b, args.min_mapped)
        if s is None:
            mapped_n = sum(1 for u in members if u in align)
            if mapped_n < args.min_mapped:
                n_skipped_coverage += 1
            else:
                n_skipped_no_edges += 1
            continue
        s["cluster_id"] = cid
        scores.append(s)

    print(f"  scored {len(scores)} clusters; skipped {n_skipped_coverage} "
          f"(< {args.min_mapped} mapped members), {n_skipped_no_edges} "
          f"(no internal edges on one side)")

    if not scores:
        print("No clusters passed filters — nothing to plot.")
        return

    for m in ("EC", "ICS", "S3"):
        vals = [s[m] for s in scores]
        print(f"  median {m}: {np.median(vals):.3f}  (mean {np.mean(vals):.3f})")

    plot(scores, "cluster_alignment.png", args.net1_name, args.net2_name, args.min_mapped)

    out = {
        "params": {
            "alignment_file": args.alignment,
            "min_mapped": args.min_mapped,
            "species_a": args.net1_name,
            "species_b": args.net2_name,
        },
        "n_clusters_total": len(clusters_a),
        "n_scored": len(scores),
        "n_skipped_coverage": n_skipped_coverage,
        "n_skipped_no_edges": n_skipped_no_edges,
        "medians": {m: float(np.median([s[m] for s in scores])) for m in ("EC", "ICS", "S3")},
        "per_cluster": scores,
    }
    with open("cluster_alignment.json", "w") as f:
        json.dump(out, f, indent=2)
    print("  wrote cluster_alignment.json")


if __name__ == "__main__":
    main()
