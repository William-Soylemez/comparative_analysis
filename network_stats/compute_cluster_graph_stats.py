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
    - average node connectivity AND average edge connectivity (the all-pairs
      local-connectivity averages that were dropped as infeasible on the full
      protein graph in compute_network_stats.py) — computed EXACTLY here
      (every C(n,2) pair of the cluster graph's largest CC, each an
      independent max-flow) by fanning the pair list out across --procs
      forked processes; this is embarrassingly parallel across pairs (see
      bench_node_connectivity.py, where this pattern was validated for the
      node case), so unlike the full protein graph, the cluster graph (2-3
      orders of magnitude smaller) is tractable this way even on a single
      multi-core node. networkx has no built-in average_edge_connectivity,
      so it's computed the same way as average_node_connectivity, just
      swapping in local_edge_connectivity + build_auxiliary_edge_connectivity
      (verified to reproduce a from-scratch, no-reuse computation exactly).

Resume: if --out already exists and has non-null avg_node_connectivity_largest_cc
and avg_edge_connectivity_largest_cc for the same --min-crossing-edges, this
is a fast no-op (nothing is rebuilt or recomputed). Pass --force to recompute
anyway.

Usage:
    python3 compute_cluster_graph_stats.py SPECIES_DIR [--min-crossing-edges N]
        [--skip STAT ...] [--out FILE] [--procs N] [--force] [--quiet]

Example:
    python3 compute_cluster_graph_stats.py ../GCF_000182965.3 --out out.json --procs 144
"""

import argparse
import glob
import json
import math
import os
import random
import sys
import time
from collections import defaultdict
from multiprocessing import get_context
from statistics import mean

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity as nx_modularity
from networkx.algorithms.connectivity import (
    build_auxiliary_edge_connectivity,
    build_auxiliary_node_connectivity,
    local_edge_connectivity,
    local_node_connectivity,
)
from networkx.algorithms.flow import build_residual_network
from tqdm import tqdm

MIN_CROSSING_EDGES = 50
ALL_SKIPPABLE = ["num_bridges", "avg_node_connectivity", "avg_edge_connectivity"]
AVG_CONNECTIVITY_SAMPLE_PAIRS = 200
# Safety net only -- not a realistic budget check anymore now that the exact
# computation is parallelized across --procs cores (see
# parallel_average_connectivity). At ~0.3-0.5 ms/pair/core this is ~30-80
# core-hours; if a cluster graph's largest CC is ever this big, something
# upstream (min-crossing-edges threshold, clustering) is probably off and
# it's worth a look before burning that much compute.
AVG_CONNECTIVITY_SAFETY_MAX_PAIRS = 20_000_000

# Module globals so forked worker processes inherit the graph + its
# auxiliary/residual structures copy-on-write, without pickling them to every
# worker (matches the approach validated in bench_node_connectivity.py).
# _CKIND picks node- vs edge-connectivity per call to parallel_average_connectivity;
# it's a module global (not a Pool argument) for the same copy-on-write-via-fork
# reason as _CG -- workers read it once at fork time via the initializer.
_CG = None
_CH = None
_CR = None
_CKIND = None


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


def _init_conn_worker():
    # Built once per forked worker, reused across every pair that worker gets.
    # _CKIND picks which flavor of auxiliary graph to build; local_node_connectivity
    # and local_edge_connectivity both accept the same auxiliary/residual reuse
    # pattern (verified to give bit-identical results to fresh per-pair calls).
    global _CH, _CR
    if _CKIND == "node":
        _CH = build_auxiliary_node_connectivity(_CG)
    else:
        _CH = build_auxiliary_edge_connectivity(_CG)
    _CR = build_residual_network(_CH, "capacity")


def _kappa(pair):
    u, v = pair
    if _CKIND == "node":
        return local_node_connectivity(_CG, u, v, auxiliary=_CH, residual=_CR)
    return local_edge_connectivity(_CG, u, v, auxiliary=_CH, residual=_CR)


def parallel_average_connectivity(g, procs, kind, quiet=False):
    """Exact average local connectivity of g (node- or edge-flavored,
    matching nx.average_node_connectivity's definition but generalized to
    edges too, since networkx has no built-in average_edge_connectivity),
    computed by fanning every C(n,2) pair out across `procs` forked processes.
    Each pair's local connectivity is an independent max-flow, so this is
    embarrassingly parallel -- see bench_node_connectivity.py, where this
    pattern was validated for the node case; the edge case reuses the same
    auxiliary/residual-network reuse trick via build_auxiliary_edge_connectivity."""
    assert kind in ("node", "edge")
    global _CG, _CKIND

    n = g.number_of_nodes()
    total_pairs = n * (n - 1) // 2
    if total_pairs == 0:
        return 0.0, {"total_pairs": 0, "procs": procs, "skipped": False, "reason": "fewer than 2 nodes"}

    if total_pairs > AVG_CONNECTIVITY_SAFETY_MAX_PAIRS:
        return None, {"total_pairs": total_pairs, "procs": procs, "skipped": True,
                       "reason": f"total_pairs exceeds safety cap {AVG_CONNECTIVITY_SAFETY_MAX_PAIRS}"}

    _CG = g
    _CKIND = kind
    nodes = list(g.nodes())

    # Quick single-core sample first, purely to report a per-pair cost and
    # projected serial-vs-parallel time -- doesn't gate anything.
    sample_n = min(AVG_CONNECTIVITY_SAMPLE_PAIRS, total_pairs)
    rng = random.Random(0)
    sample_pairs = set()
    while len(sample_pairs) < sample_n:
        u, v = rng.sample(nodes, 2)
        sample_pairs.add((u, v) if u < v else (v, u))
    _init_conn_worker()
    t0 = time.time()
    for pr in sample_pairs:
        _kappa(pr)
    ms_per_pair = (time.time() - t0) / sample_n * 1000

    if not quiet:
        print(f"  avg_{kind}_connectivity: {ms_per_pair:.2f} ms/pair (1 core), "
              f"{total_pairs:,} pairs total -> ~{human(total_pairs*ms_per_pair/1000):s} serial, "
              f"~{human(total_pairs*ms_per_pair/1000/procs):s} on {procs} cores", file=sys.stderr)

    all_pairs = [(nodes[i], nodes[j]) for i in range(n) for j in range(i + 1, n)]

    t0 = time.time()
    ctx = get_context("fork")
    with ctx.Pool(processes=procs, initializer=_init_conn_worker) as pool:
        kappas = pool.map(_kappa, all_pairs, chunksize=max(1, len(all_pairs) // (procs * 4)))
    wall = time.time() - t0

    if not quiet:
        print(f"  computed exactly: {total_pairs:,} pairs in {human(wall)} on {procs} cores", file=sys.stderr)

    meta = {
        "sample_pairs": sample_n, "ms_per_pair": ms_per_pair, "total_pairs": total_pairs,
        "procs": procs, "parallel_seconds": wall, "skipped": False,
    }
    return mean(kappas), meta


def human(seconds):
    if seconds < 90:
        return f"{seconds:.1f}s"
    if seconds < 5400:
        return f"{seconds/60:.1f}min"
    return f"{seconds/3600:.1f}h"


# Fields that require the cluster graph to be rebuilt but are otherwise
# cheap (cluster graphs are 2-3 orders of magnitude smaller than the protein
# graph). Tracked separately from avg_node_connectivity_largest_cc so a
# species that already has a valid (expensive, exact) connectivity value
# never has it recomputed just because a new cheap field (like modularity)
# was added later -- only the new field gets backfilled.
CHEAP_FIELDS = [
    "node_count", "edge_count", "num_connected_components", "largest_cc_size",
    "global_clustering_coefficient", "num_bridges", "node_connectivity_largest_cc",
    "edge_connectivity_largest_cc", "diameter_largest_cc", "modularity",
]


def compute_stats(species_dir, min_crossing, skip, quiet=False, procs=None, existing=None):
    """existing: a previously-computed result dict for this species (same
    schema). Lets this be safely rerun to backfill just what's missing --
    the expensive avg_node_connectivity is never recomputed once it has a
    real value, and cheap fields (e.g. added in a later version of this
    script) are filled in without re-touching avg_node_connectivity."""
    acc, network_path, clusters_path = find_species_files(species_dir)
    existing = dict(existing) if existing else None
    same_threshold = bool(existing) and existing.get("min_crossing_edges") == min_crossing

    have_node_conn = same_threshold and existing.get("avg_node_connectivity_largest_cc") is not None
    have_edge_conn = same_threshold and existing.get("avg_edge_connectivity_largest_cc") is not None
    have_cheap = same_threshold and all(k in existing for k in CHEAP_FIELDS)

    if have_node_conn and have_edge_conn and have_cheap:
        if not quiet:
            print(f"species: {acc} -- already complete, nothing to do", file=sys.stderr)
        return existing

    if not quiet:
        print(f"species: {acc}", file=sys.stderr)
    procs = procs or os.cpu_count()

    p2c, all_cluster_ids = load_protein_to_clusters(clusters_path)
    g = build_cluster_graph(network_path, p2c, all_cluster_ids, min_crossing, quiet=quiet)

    components = list(nx.connected_components(g))
    largest_cc = max(components, key=len) if components else set()
    sub = g.subgraph(largest_cc)

    result = dict(existing) if same_threshold and existing else {"species": acc, "min_crossing_edges": min_crossing}

    if not have_cheap:
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

    if "avg_node_connectivity" in skip:
        result["avg_node_connectivity_largest_cc"] = None
        result["avg_node_connectivity_meta"] = {"skipped": True, "reason": "--skip"}
    elif not have_node_conn:
        val, meta = parallel_average_connectivity(sub, procs, "node", quiet=quiet)
        result["avg_node_connectivity_largest_cc"] = val
        result["avg_node_connectivity_meta"] = meta
    # else: have_node_conn already True, its value is preserved via the dict(existing) copy above.

    if "avg_edge_connectivity" in skip:
        result["avg_edge_connectivity_largest_cc"] = None
        result["avg_edge_connectivity_meta"] = {"skipped": True, "reason": "--skip"}
    elif not have_edge_conn:
        val, meta = parallel_average_connectivity(sub, procs, "edge", quiet=quiet)
        result["avg_edge_connectivity_largest_cc"] = val
        result["avg_edge_connectivity_meta"] = meta
    # else: have_edge_conn already True, its value is preserved via the dict(existing) copy above.

    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../GCF_000182965.3)")
    ap.add_argument("--min-crossing-edges", type=int, default=MIN_CROSSING_EDGES,
                     help=f"Minimum crossing protein-protein edges for a cluster-graph edge (default {MIN_CROSSING_EDGES})")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output): {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output JSON path (default: print to stdout)")
    ap.add_argument("--procs", type=int, default=os.cpu_count(),
                     help="worker processes for avg_node_connectivity (default: all cores)")
    ap.add_argument("--force", action="store_true",
                     help="Recompute avg_node_connectivity even if --out already has a value")
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
                            quiet=args.quiet, procs=args.procs, existing=existing)

    out_json = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_json + "\n")
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        print(out_json)


if __name__ == "__main__":
    main()
