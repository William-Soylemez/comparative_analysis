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
    out_file_abs = os.path.abspath(out_file)

    results = []
    for path in sorted(glob.glob(f"{stats_dir}/*_stats.json")):
        # The default out_file name (all_species_stats.json) itself matches
        # this *_stats.json glob -- without this check, any rerun (all three
        # batch scripts call this as their last step) would re-ingest the
        # PREVIOUS combined output as if it were one more species record,
        # nesting it one level deeper each time (this bit us: 6 levels deep
        # after repeated reruns). Comparing absolute paths so this works
        # regardless of what out_file is actually named or where it lives.
        if os.path.abspath(path) == out_file_abs:
            continue
        with open(path) as f:
            loaded = json.load(f)
        if not isinstance(loaded, dict):
            print(f"skipping {path}: expected a single species object, got {type(loaded).__name__} "
                  f"(this file may itself be corrupted -- e.g. leftover nested aggregate output)",
                  file=sys.stderr)
            continue
        results.append(loaded)

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
