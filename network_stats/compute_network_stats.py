#!/usr/bin/env python3
"""
Global network statistics for a single PHILHARMONIC species run.

Reads `<acc>_network.positive.tsv` (and `<acc>_clusters.json` for the cluster
count) from a species directory and reports, as JSON:

    - density
    - node count
    - edge count
    - number of connected components
    - number of nodes outside the largest connected component
    - median degree
    - number of PHILHARMONIC clusters
    - number of bridges
    - global clustering coefficient (transitivity)

Note: average node/edge connectivity (networkx's `average_node_connectivity`
and the edge-connectivity analog) were evaluated and dropped. They require
local connectivity between every pair of nodes via max-flow — benchmarked at
~425ms/pair, i.e. ~51 days for the smallest species here (4,555 nodes) even
restricted to the largest connected component. Not viable at this scale.

Single-threaded by design: this is the per-species unit of work for a later
SLURM array that runs many species in parallel.

Usage:
    python3 compute_network_stats.py SPECIES_DIR [--skip STAT ...] [--out FILE]

Example:
    python3 compute_network_stats.py ../GCF_000182965.3 --skip num_bridges
"""

import argparse
import glob
import json
import os
import sys
from statistics import median

import networkx as nx
from tqdm import tqdm

ALL_SKIPPABLE = ["num_bridges"]


def find_species_files(species_dir):
    species_dir = species_dir.rstrip("/")
    net_hits = glob.glob(f"{species_dir}/*_network.positive.tsv")
    if not net_hits:
        sys.exit(f"No *_network.positive.tsv found in {species_dir}")
    network_path = net_hits[0]
    acc = os.path.basename(network_path).replace("_network.positive.tsv", "")

    cluster_hits = glob.glob(f"{species_dir}/*_clusters.json")
    clusters_path = cluster_hits[0] if cluster_hits else None
    return acc, network_path, clusters_path


def load_graph(network_path):
    """Read an undirected, deduplicated, self-loop-free edge list into a Graph."""
    g = nx.Graph()
    with open(network_path) as f:
        lines = f.readlines()
    for line in tqdm(lines, desc="loading edges", unit="edge"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        a, b = parts[0], parts[1]
        if a == b:
            continue
        g.add_edge(a, b)
    return g


def count_clusters(clusters_path):
    if clusters_path is None:
        return None
    with open(clusters_path) as f:
        return len(json.load(f))


def count_bridges(g, components):
    """Sum bridges within each connected component (nx.bridges requires connectivity)."""
    total = 0
    for comp in tqdm(components, desc="bridges per component", unit="component"):
        sub = g.subgraph(comp)
        if sub.number_of_edges() == 0:
            continue
        total += sum(1 for _ in nx.bridges(sub))
    return total


def compute_stats(species_dir, skip):
    acc, network_path, clusters_path = find_species_files(species_dir)
    print(f"species: {acc}", file=sys.stderr)

    g = load_graph(network_path)
    degrees = [d for _, d in g.degree()]
    components = list(nx.connected_components(g))
    largest_cc = max(components, key=len)

    result = {
        "species": acc,
        "node_count": g.number_of_nodes(),
        "edge_count": g.number_of_edges(),
        "density": nx.density(g),
        "num_connected_components": len(components),
        "nodes_outside_largest_cc": g.number_of_nodes() - len(largest_cc),
        "median_degree": median(degrees) if degrees else 0,
        "num_philharmonic_clusters": count_clusters(clusters_path),
        "global_clustering_coefficient": nx.transitivity(g),
        "largest_cc_size": len(largest_cc),
    }

    if "num_bridges" in skip:
        result["num_bridges"] = None
    else:
        result["num_bridges"] = count_bridges(g, components)

    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../GCF_000182965.3)")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output) to save time: {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output JSON path (default: print to stdout)")
    args = ap.parse_args()

    result = compute_stats(args.species_dir, set(args.skip))

    out_json = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_json + "\n")
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(out_json)


if __name__ == "__main__":
    main()
