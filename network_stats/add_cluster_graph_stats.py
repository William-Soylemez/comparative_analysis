#!/usr/bin/env python3
"""
Computes the cluster-graph metrics from compute_cluster_graph_stats.py and
merges them into an existing per-species stats JSON (as written by
compute_network_stats.py), nested under a "cluster_graph" key.

Nesting matters here, not just style: compute_network_stats.py's protein-
graph stats and compute_cluster_graph_stats.py's cluster-graph stats use the
same field names (e.g. "node_count", "edge_count", "global_clustering_
coefficient") for two different graphs. A flat merge would silently
overwrite one graph's numbers with the other's.

Usage:
    python3 add_cluster_graph_stats.py SPECIES_DIR STATS_JSON
        [--min-crossing-edges N] [--skip STAT ...] [--out FILE] [--quiet]

With no --out, STATS_JSON is updated in place.

Example:
    python3 add_cluster_graph_stats.py ../GCF_000182965.3 network_stats_out/GCF_000182965.3_stats.json
"""

import argparse
import json
import os
import sys

from compute_cluster_graph_stats import compute_stats, MIN_CROSSING_EDGES, ALL_SKIPPABLE


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir", help="Path to a species directory (e.g. ../GCF_000182965.3)")
    ap.add_argument("stats_json", help="Existing per-species stats JSON to update (from compute_network_stats.py)")
    ap.add_argument("--min-crossing-edges", type=int, default=MIN_CROSSING_EDGES,
                     help=f"Minimum crossing protein-protein edges for a cluster-graph edge (default {MIN_CROSSING_EDGES})")
    ap.add_argument("--skip", nargs="+", default=[], choices=ALL_SKIPPABLE,
                     help=f"Stats to skip (set to null in output): {ALL_SKIPPABLE}")
    ap.add_argument("--out", default=None, help="Output path (default: overwrite stats_json in place)")
    ap.add_argument("--quiet", action="store_true", help="Disable progress bars")
    args = ap.parse_args()

    if not os.path.exists(args.stats_json):
        sys.exit(f"{args.stats_json} does not exist -- run compute_network_stats.py first")

    with open(args.stats_json) as f:
        existing = json.load(f)

    cluster_stats = compute_stats(args.species_dir, args.min_crossing_edges, set(args.skip), quiet=args.quiet)
    cluster_stats.pop("species", None)  # already present at the top level of the stats JSON

    existing["cluster_graph"] = cluster_stats

    out_path = args.out or args.stats_json
    with open(out_path, "w") as f:
        json.dump(existing, f, indent=2)
        f.write("\n")
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
