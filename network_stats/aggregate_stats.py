#!/usr/bin/env python3
"""
Combine per-species network-stats JSON files (as written by
compute_network_stats.py --out) into one JSON list.

Usage:
    python3 aggregate_stats.py STATS_DIR [OUT_FILE]
"""

import glob
import json
import os
import sys


def main():
    stats_dir = sys.argv[1].rstrip("/")
    out_file = sys.argv[2] if len(sys.argv) > 2 else f"{stats_dir}/all_species_stats.json"

    results = []
    for path in sorted(glob.glob(f"{stats_dir}/*_stats.json")):
        with open(path) as f:
            results.append(json.load(f))

    # Atomic write: run_species_batch.sh, run_node_connectivity.sh, and
    # run_edge_connectivity.sh each call this as their last step, and can
    # legitimately be running at the same time against the same OUT_DIR (see
    # their headers) -- writing straight to out_file with plain open(..., "w")
    # would let two concurrent calls interleave/truncate each other's output
    # mid-write, corrupting the file. Writing to a temp file in the same
    # directory and os.replace()-ing it into place means whichever call
    # finishes last atomically wins with a fully-valid (if possibly slightly
    # stale relative to the other, still-running call) JSON file -- never a
    # half-written or interleaved one. Run this once more manually after all
    # jobs you care about have finished for a guaranteed-fresh combined file.
    tmp_file = f"{out_file}.tmp.{os.getpid()}"
    with open(tmp_file, "w") as f:
        json.dump(results, f, indent=2)
    os.replace(tmp_file, out_file)
    print(f"wrote {out_file} ({len(results)} species)")


if __name__ == "__main__":
    main()
