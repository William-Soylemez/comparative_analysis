#!/usr/bin/env python3
"""
Combine per-species network-stats JSON files (as written by
compute_network_stats.py --out) into one JSON list.

Usage:
    python3 aggregate_stats.py STATS_DIR [OUT_FILE]
"""

import glob
import json
import sys


def main():
    stats_dir = sys.argv[1].rstrip("/")
    out_file = sys.argv[2] if len(sys.argv) > 2 else f"{stats_dir}/all_species_stats.json"

    results = []
    for path in sorted(glob.glob(f"{stats_dir}/*_stats.json")):
        with open(path) as f:
            results.append(json.load(f))

    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"wrote {out_file} ({len(results)} species)")


if __name__ == "__main__":
    main()
