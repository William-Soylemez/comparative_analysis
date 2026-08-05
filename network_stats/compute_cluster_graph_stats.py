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
    - diameter (exact, of the largest connected component) — unlike
      compute_network_stats.py's approximate protein-graph diameter, the
      cluster graph is small enough for networkx's exact nx.diameter
    - modularity of the cluster graph's own community structure, via
      networkx's greedy-modularity communities (there's no pre-existing
      "cluster of clusters" partition to score the way compute_network_stats.py
      scores the protein graph against its PHILHARMONIC clusters, so this
      instead asks how much community structure the cluster graph has on
      its own terms)

REMOVED (for now): average node/edge connectivity — the all-pairs
local-connectivity averages over the cluster graph's largest CC. These were
computed here via a parallel fork-pool over every C(n,2) pair, orchestrated by
run_node_connectivity.sh / run_edge_connectivity.sh. That whole scheme was
janky (two independent jobs racing on one output file; see the git history and
the README) and has been pulled out to be redesigned from scratch. Everything
this script now computes is cheap and single-threaded.

Resume: if --out already exists with all fields for the same
--min-crossing-edges, this is a fast no-op. Pass --force to recompute anyway.

Usage:
    python3 compute_cluster_graph_stats.py SPECIES_DIR [--min-crossing-edges N]
        [--skip STAT ...] [--out FILE] [--force] [--quiet]

Example:
    python3 compute_cluster_graph_stats.py ../input/GCF_000182965.3 --out out.json
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity as nx_modularity
from tqdm import tqdm

MIN_CROSSING_EDGES = 50
ALL_SKIPPABLE = ["num_bridges"]


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


def compute_cluster_diameter(sub, quiet=False):
    """Exact diameter of the largest CC -- unlike the protein graph, cluster
    graphs here (hundreds to ~1-2k nodes) are small enough that nx.diameter
    (all-pairs BFS) is cheap; no need for the approximation used for the
    much larger protein graph in compute_network_stats.py."""
    if sub.number_of_nodes() < 2:
        return 0
    if not quiet:
        print("  computing cluster graph diameter (exact) ...", file=sys.stderr)
    return nx.diameter(sub)


def compute_cluster_modularity(g, quiet=False):
    """Modularity of the cluster graph's own greedy-modularity communities.
    There's no pre-existing partition to score it against (unlike the protein
    graph, which is scored against the actual PHILHARMONIC clusters), so this
    asks how much community structure the cluster graph has on its own terms."""
    if g.number_of_edges() == 0:
        return None
    if not quiet:
        print("  computing cluster graph modularity (greedy communities) ...", file=sys.stderr)
    communities = greedy_modularity_communities(g)
    return nx_modularity(g, communities)


# Every field this script computes -- used for the field-level resume check.
# All are cheap and single-threaded (cluster graphs are 2-3 orders of
# magnitude smaller than the protein graph).
CHEAP_FIELDS = [
    "node_count", "edge_count", "num_connected_components", "largest_cc_size",
    "global_clustering_coefficient", "num_bridges", "node_connectivity_largest_cc",
    "edge_connectivity_largest_cc", "diameter_largest_cc", "modularity",
]


def compute_stats(species_dir, min_crossing, skip, quiet=False, existing=None):
    """existing: a previously-computed result dict for this species (same
    schema). Lets this be rerun to backfill just what's missing -- a field
    added in a later version gets filled in without recomputing a complete
    file."""
    acc, network_path, clusters_path = find_species_files(species_dir)
    existing = dict(existing) if existing else None
    same_threshold = bool(existing) and existing.get("min_crossing_edges") == min_crossing

    if same_threshold and all(k in existing for k in CHEAP_FIELDS):
        if not quiet:
            print(f"species: {acc} -- already complete, nothing to do", file=sys.stderr)
        return existing

    if not quiet:
        print(f"species: {acc}", file=sys.stderr)

    p2c, all_cluster_ids = load_protein_to_clusters(clusters_path)
    g = build_cluster_graph(network_path, p2c, all_cluster_ids, min_crossing, quiet=quiet)

    components = list(nx.connected_components(g))
    largest_cc = max(components, key=len) if components else set()
    sub = g.subgraph(largest_cc)

    result = dict(existing) if same_threshold and existing else {"species": acc, "min_crossing_edges": min_crossing}
    result["node_count"] = g.number_of_nodes()
    result["edge_count"] = g.number_of_edges()
    result["num_connected_components"] = len(components)
    result["largest_cc_size"] = len(largest_cc)
    result["global_clustering_coefficient"] = nx.transitivity(g)
    result["num_bridges"] = None if "num_bridges" in skip else count_bridges(g, components, quiet=quiet)
    result["node_connectivity_largest_cc"] = nx.node_connectivity(sub) if sub.number_of_nodes() > 1 else None
    result["edge_connectivity_largest_cc"] = nx.edge_connectivity(sub) if sub.number_of_nodes() > 1 else None
    result["diameter_largest_cc"] = compute_cluster_diameter(sub, quiet=quiet)
    result["modularity"] = compute_cluster_modularity(g, quiet=quiet)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../input/GCF_000182965.3)")
    ap.add_argument("--min-crossing-edges", type=int, default=MIN_CROSSING_EDGES,
                     help=f"Minimum crossing protein-protein edges for a cluster-graph edge (default {MIN_CROSSING_EDGES})")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output): {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output JSON path (default: print to stdout)")
    ap.add_argument("--force", action="store_true",
                     help="Recompute even if --out already has complete stats")
    ap.add_argument("--quiet", action="store_true", help="Disable progress bars")
    args = ap.parse_args()

    existing = None
    if args.out and not args.force and os.path.exists(args.out):
        with open(args.out) as f:
            try:
                existing = json.load(f)
            except json.JSONDecodeError:
                existing = None

    result = compute_stats(args.species_dir, args.min_crossing_edges, set(args.skip),
                            quiet=args.quiet, existing=existing)

    out_json = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_json + "\n")
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(out_json)


if __name__ == "__main__":
    main()
