#!/usr/bin/env python3
"""
Cluster-graph statistics for a single PHILHARMONIC species run — a "cluster
of clusters" graph, where each node is a PHILHARMONIC cluster (not a
protein), and an edge exists between two clusters if at least
MIN_CROSSING_EDGES protein-protein edges cross between their members in the
underlying PPI network (`<acc>_network.positive.tsv`).

Built from scratch here rather than reusing the PHILHARMONIC-provided
`<acc>_cluster_graph.tsv` — spot-checking one species (GCF_000182965.3)
found its weight for a specific cluster pair (65) didn't match the actual
crossing-edge count between those clusters' full `clusters.json` member
lists (133), most likely because it was built from core (pre-ReCIPE) rather
than full cluster membership. Rebuilding directly from `clusters.json` +
the network file keeps this consistent with every other analysis in this
project, which all use `clusters.json`'s full `members` list.

Proteins that belong to more than one cluster (ReCIPE can re-add a hub to
several clusters) are NOT collapsed to a single cluster — an edge between
proteins u and v counts as crossing for every (cluster of u, cluster of v)
combination.

Reports, as JSON:
    - node count (= number of PHILHARMONIC clusters), edge count of the
      cluster graph
    - number of connected components, largest CC size
    - global clustering coefficient (transitivity)
    - number of bridges
    - node connectivity, edge connectivity (standard graph invariants —
      minimum nodes/edges whose removal disconnects the graph), computed on
      the largest connected component
    - average node connectivity (the all-pairs local-connectivity average
      that was dropped as infeasible on the full protein graph in
      compute_network_stats.py) — attempted here since the cluster graph is
      2-3 orders of magnitude smaller; benchmarked on a small sample first
      and skipped with a note if the extrapolated cost is still too high

Usage:
    python3 compute_cluster_graph_stats.py SPECIES_DIR [--min-crossing-edges N]
        [--skip STAT ...] [--out FILE] [--quiet]

Example:
    python3 compute_cluster_graph_stats.py ../GCF_000182965.3
"""

import argparse
import glob
import json
import os
import sys
import time
from collections import defaultdict

import networkx as nx
from tqdm import tqdm

MIN_CROSSING_EDGES = 50
ALL_SKIPPABLE = ["num_bridges", "avg_node_connectivity"]
AVG_CONNECTIVITY_SAMPLE_PAIRS = 200
AVG_CONNECTIVITY_BUDGET_SECONDS = 1800  # skip if extrapolated full cost exceeds this (30 min)


def find_species_files(species_dir):
    species_dir = species_dir.rstrip("/")
    net_hits = glob.glob(f"{species_dir}/*_network.positive.tsv")
    if not net_hits:
        sys.exit(f"No *_network.positive.tsv found in {species_dir}")
    network_path = net_hits[0]
    acc = os.path.basename(network_path).replace("_network.positive.tsv", "")

    cluster_hits = glob.glob(f"{species_dir}/*_clusters.json")
    if not cluster_hits:
        sys.exit(f"No *_clusters.json found in {species_dir}")
    return acc, network_path, cluster_hits[0]


def load_protein_to_clusters(clusters_path):
    """protein_id -> set of cluster ids it's a member of (usually size 1)."""
    with open(clusters_path) as f:
        clusters = json.load(f)
    p2c = defaultdict(set)
    for cid, c in clusters.items():
        for m in c["members"]:
            p2c[m].add(cid)
    return p2c, set(clusters.keys())


def build_cluster_graph(network_path, p2c, all_cluster_ids, min_crossing, quiet=False):
    crossing = defaultdict(int)
    with open(network_path) as f:
        lines = f.readlines()
    for line in tqdm(lines, desc="counting crossing edges", unit="edge", disable=quiet):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        u, v = parts[0], parts[1]
        if u == v:
            continue
        cu, cv = p2c.get(u), p2c.get(v)
        if not cu or not cv:
            continue
        for a in cu:
            for b in cv:
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                crossing[key] += 1

    g = nx.Graph()
    g.add_nodes_from(all_cluster_ids)
    for (a, b), count in crossing.items():
        if count >= min_crossing:
            g.add_edge(a, b, weight=count)
    return g


def count_bridges(g, components, quiet=False):
    total = 0
    for comp in tqdm(components, desc="bridges per component", unit="component", disable=quiet):
        sub = g.subgraph(comp)
        if sub.number_of_edges() == 0:
            continue
        total += sum(1 for _ in nx.bridges(sub))
    return total


def estimate_and_maybe_compute_avg_connectivity(g, quiet=False):
    """Benchmark local_node_connectivity on a sample of pairs; only run the
    full average_node_connectivity if the extrapolated cost fits the budget."""
    import random
    from networkx.algorithms.connectivity import local_node_connectivity

    n = g.number_of_nodes()
    total_pairs = n * (n - 1) // 2
    if total_pairs == 0:
        return 0.0, {"sampled": False, "reason": "fewer than 2 nodes"}

    nodes = list(g.nodes())
    sample_n = min(AVG_CONNECTIVITY_SAMPLE_PAIRS, total_pairs)
    rng = random.Random(0)
    sample_pairs = set()
    while len(sample_pairs) < sample_n:
        u, v = rng.sample(nodes, 2)
        key = (u, v) if u < v else (v, u)
        sample_pairs.add(key)

    t0 = time.time()
    for u, v in sample_pairs:
        local_node_connectivity(g, u, v)
    elapsed = time.time() - t0
    per_pair = elapsed / sample_n
    extrapolated = per_pair * total_pairs

    if not quiet:
        print(f"  avg_node_connectivity: {per_pair*1000:.2f}ms/pair sampled, "
              f"extrapolated {extrapolated:.1f}s for all {total_pairs} pairs "
              f"(budget {AVG_CONNECTIVITY_BUDGET_SECONDS}s)", file=sys.stderr)

    meta = {"sampled": True, "sample_pairs": sample_n, "ms_per_pair": per_pair * 1000,
            "extrapolated_seconds": extrapolated, "total_pairs": total_pairs}

    if extrapolated > AVG_CONNECTIVITY_BUDGET_SECONDS:
        meta["skipped"] = True
        return None, meta

    if not quiet:
        print("  computing average_node_connectivity ...", file=sys.stderr)
    meta["skipped"] = False
    return nx.average_node_connectivity(g), meta


def compute_stats(species_dir, min_crossing, skip, quiet=False):
    acc, network_path, clusters_path = find_species_files(species_dir)
    if not quiet:
        print(f"species: {acc}", file=sys.stderr)

    p2c, all_cluster_ids = load_protein_to_clusters(clusters_path)
    g = build_cluster_graph(network_path, p2c, all_cluster_ids, min_crossing, quiet=quiet)

    components = list(nx.connected_components(g))
    largest_cc = max(components, key=len) if components else set()
    sub = g.subgraph(largest_cc)

    result = {
        "species": acc,
        "min_crossing_edges": min_crossing,
        "node_count": g.number_of_nodes(),
        "edge_count": g.number_of_edges(),
        "num_connected_components": len(components),
        "largest_cc_size": len(largest_cc),
        "global_clustering_coefficient": nx.transitivity(g),
    }

    if "num_bridges" in skip:
        result["num_bridges"] = None
    else:
        result["num_bridges"] = count_bridges(g, components, quiet=quiet)

    result["node_connectivity_largest_cc"] = nx.node_connectivity(sub) if sub.number_of_nodes() > 1 else None
    result["edge_connectivity_largest_cc"] = nx.edge_connectivity(sub) if sub.number_of_nodes() > 1 else None

    if "avg_node_connectivity" in skip:
        result["avg_node_connectivity_largest_cc"] = None
        result["avg_node_connectivity_meta"] = {"skipped": True, "reason": "--skip"}
    else:
        val, meta = estimate_and_maybe_compute_avg_connectivity(sub, quiet=quiet)
        result["avg_node_connectivity_largest_cc"] = val
        result["avg_node_connectivity_meta"] = meta

    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../GCF_000182965.3)")
    ap.add_argument("--min-crossing-edges", type=int, default=MIN_CROSSING_EDGES,
                     help=f"Minimum crossing protein-protein edges for a cluster-graph edge (default {MIN_CROSSING_EDGES})")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output): {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output JSON path (default: print to stdout)")
    ap.add_argument("--quiet", action="store_true", help="Disable progress bars")
    args = ap.parse_args()

    result = compute_stats(args.species_dir, args.min_crossing_edges, set(args.skip), quiet=args.quiet)

    out_json = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_json + "\n")
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(out_json)


if __name__ == "__main__":
    main()
