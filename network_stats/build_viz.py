#!/usr/bin/env python3
"""
Build the species-stats visualization by embedding the aggregated stats JSON
into the HTML template.

Reads results/all_species_data.json (the output of aggregate_stats.py), splices
it into species_viz.template.html at the `__EMBEDDED_DATA__` marker, and writes
results/species_viz.html.

The embedded copy is what makes the page work on a plain double-click
(file://). When the page is instead *served* (e.g. `python3 -m http.server`),
it also fetches all_species_data.json from the same directory at load time and
uses that live copy if present -- so the embedded data is just the offline
fallback, and re-running this script keeps that fallback in sync with the
current results.

Usage (from this directory):
    ../venv/bin/python build_viz.py
    ../venv/bin/python build_viz.py --data results/all_species_data.json --out results/species_viz.html
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MARKER = "__EMBEDDED_DATA__"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.join(HERE, "results", "all_species_data.json"),
                    help="aggregated stats JSON to embed (default results/all_species_data.json)")
    ap.add_argument("--template", default=os.path.join(HERE, "species_viz.template.html"),
                    help="HTML template with the %s marker" % MARKER)
    ap.add_argument("--out", default=os.path.join(HERE, "results", "species_viz.html"),
                    help="output HTML path (default results/species_viz.html)")
    args = ap.parse_args()

    with open(args.data) as f:
        data = json.load(f)
    if not isinstance(data, list):
        sys.exit(f"{args.data}: expected a JSON array of species objects, got {type(data).__name__}")

    template = open(args.template).read()
    if MARKER not in template:
        sys.exit(f"{args.template}: missing {MARKER} marker")

    # json.dumps output is a valid JS array literal; drop it straight in.
    html = template.replace(MARKER, json.dumps(data))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(html)
    print(f"wrote {args.out} ({len(data)} species embedded)")


if __name__ == "__main__":
    main()
