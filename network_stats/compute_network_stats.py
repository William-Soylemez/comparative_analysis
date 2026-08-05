#!/usr/bin/env python3
"""
Global network statistics for a single PHILHARMONIC species run.

Reads `<acc>_network.positive.tsv` (and `<acc>_clusters.json`) from a
species directory and reports, as JSON:

    - density, node count, edge count
    - number of connected components, nodes outside the largest one
    - median degree
    - number of PHILHARMONIC clusters
    - number of bridges
    - global clustering coefficient (transitivity)
    - cluster_graph: a nested block of stats on the "cluster of clusters"
      graph (see compute_cluster_graph_stats.py) — clustering coefficient,
      bridges, node/edge connectivity, diameter, and modularity, all on the
      cluster-level graph, not the protein-level one. Nested rather than
      flattened because both graphs use the same field names (node_count,
      edge_count, ...) for different things.

Field-level resume: if --out already exists and parses as JSON, only the
keys MISSING from it are (re)computed — an already-complete file is a fast
no-op (just a key check, no network file read). This means
compute_network_stats.py can just be rerun for everything after adding a
new stat (like cluster_graph, added after the basic protein-graph stats
already existed for some species) without needing a separate script or
deleting existing output first. Pass --force to ignore existing output and
recompute everything.

Note: diameter and modularity are computed only for the cluster graph (see
compute_cluster_graph_stats.py), not the protein graph.

Note: average node/edge connectivity (all-pairs local connectivity) is NOT
computed. On the protein graph it's infeasible (~425ms/pair, ~51 days for the
smallest species here, 4,555 nodes). A parallelized cluster-graph version used
to exist, but its orchestration was janky and it has been removed pending a
redesign — see the README.

Single-threaded by design: this is the per-species unit of work for a later
SLURM array that runs many species in parallel.

Usage:
    python3 compute_network_stats.py SPECIES_DIR [--skip STAT ...] [--out FILE] [--force]

Example:
    python3 compute_network_stats.py ../input/GCF_000182965.3 --out out/GCF_000182965.3_stats.json
    # rerun later after a new stat is added -- only computes what's missing:
    python3 compute_network_stats.py ../input/GCF_000182965.3 --out out/GCF_000182965.3_stats.json
"""

import argparse
import glob
import json
import os
import sys
from statistics import median

import networkx as nx
from tqdm import tqdm

from compute_cluster_graph_stats import (
    compute_stats as compute_cluster_graph_stats,
    MIN_CROSSING_EDGES as CLUSTER_MIN_CROSSING_EDGES,
    CHEAP_FIELDS as CLUSTER_CHEAP_FIELDS,
)

ALL_SKIPPABLE = ["num_bridges", "cluster_graph", "cluster_num_bridges"]

EXPECTED_KEYS = [
    "node_count", "edge_count", "density", "num_connected_components",
    "nodes_outside_largest_cc", "median_degree", "num_philharmonic_clusters",
    "global_clustering_coefficient", "largest_cc_size", "num_bridges",
    "cluster_graph",
]


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


def load_graph(network_path, quiet=False):
    """Read an undirected, deduplicated, self-loop-free edge list into a Graph."""
    g = nx.Graph()
    with open(network_path) as f:
        lines = f.readlines()
    for line in tqdm(lines, desc="loading edges", unit="edge", disable=quiet):
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


def count_bridges(g, components, quiet=False):
    """Sum bridges within each connected component (nx.bridges requires connectivity)."""
    total = 0
    for comp in tqdm(components, desc="bridges per component", unit="component", disable=quiet):
        sub = g.subgraph(comp)
        if sub.number_of_edges() == 0:
            continue
        total += sum(1 for _ in nx.bridges(sub))
    return total


def compute_stats(species_dir, skip, quiet=False, existing=None):
    existing = dict(existing) if existing else {}
    acc, network_path, clusters_path = find_species_files(species_dir)

    missing = [k for k in EXPECTED_KEYS if k not in existing]

    # cluster_graph is present but missing a field added in a later version of
    # compute_cluster_graph_stats.py (e.g. diameter_largest_cc/modularity) --
    # retry it to backfill just that, unless THIS run is skipping cluster_graph.
    # Passing the existing cluster_graph dict through as `existing` lets
    # compute_cluster_graph_stats fill only what's missing.
    cg = existing.get("cluster_graph")
    if "cluster_graph" not in missing and cg is not None and "cluster_graph" not in skip:
        if any(k not in cg for k in CLUSTER_CHEAP_FIELDS):
            missing = missing + ["cluster_graph"]

    if not missing:
        if not quiet:
            print(f"species: {acc} -- already complete, nothing to do", file=sys.stderr)
        return existing

    if not quiet:
        print(f"species: {acc} -- computing missing: {missing}", file=sys.stderr)

    result = dict(existing)
    result["species"] = acc

    protein_level_missing = [k for k in missing if k != "cluster_graph"]
    if protein_level_missing:
        g = load_graph(network_path, quiet=quiet)
        degrees = [d for _, d in g.degree()]
        components = list(nx.connected_components(g))
        largest_cc = max(components, key=len)

        if "node_count" in missing:
            result["node_count"] = g.number_of_nodes()
        if "edge_count" in missing:
            result["edge_count"] = g.number_of_edges()
        if "density" in missing:
            result["density"] = nx.density(g)
        if "num_connected_components" in missing:
            result["num_connected_components"] = len(components)
        if "nodes_outside_largest_cc" in missing:
            result["nodes_outside_largest_cc"] = g.number_of_nodes() - len(largest_cc)
        if "median_degree" in missing:
            result["median_degree"] = median(degrees) if degrees else 0
        if "num_philharmonic_clusters" in missing:
            result["num_philharmonic_clusters"] = count_clusters(clusters_path)
        if "global_clustering_coefficient" in missing:
            result["global_clustering_coefficient"] = nx.transitivity(g)
        if "largest_cc_size" in missing:
            result["largest_cc_size"] = len(largest_cc)

        if "num_bridges" in missing:
            result["num_bridges"] = None if "num_bridges" in skip else count_bridges(g, components, quiet=quiet)

    if "cluster_graph" in missing:
        if "cluster_graph" in skip:
            result["cluster_graph"] = None
        else:
            cluster_skip = set()
            if "cluster_num_bridges" in skip:
                cluster_skip.add("num_bridges")
            cs = compute_cluster_graph_stats(species_dir, CLUSTER_MIN_CROSSING_EDGES, cluster_skip,
                                              quiet=quiet, existing=cg)
            cs.pop("species", None)
            result["cluster_graph"] = cs

    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../input/GCF_000182965.3)")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output) to save time: {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output JSON path (default: print to stdout)")
    ap.add_argument("--force", action="store_true",
                     help="Recompute everything even if --out already has complete stats")
    ap.add_argument("--quiet", action="store_true",
                     help="Disable progress bars (use when running many species in parallel)")
    args = ap.parse_args()

    existing = {}
    if args.out and not args.force and os.path.exists(args.out):
        try:
            with open(args.out) as f:
                existing = json.load(f)
        except (json.JSONDecodeError, OSError):
            existing = {}

    result = compute_stats(args.species_dir, set(args.skip), quiet=args.quiet, existing=existing)

    if args.out:
        with open(args.out, "w") as f:
            f.write(json.dumps(result, indent=2) + "\n")
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
