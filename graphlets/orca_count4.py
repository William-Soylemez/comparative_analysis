#!/usr/bin/env python3
"""
Count connected n=4 graphlets for each species' PPI network using ORCA
(fast combinatorial orbit counting), which gives the SAME exact counts as
exhaustive enumeration without walking every 4-subgraph.

Pipeline per species:
  1. read <acc>_network.positive.tsv, remap protein ids -> 0..n-1,
     dedupe undirected edges, drop self-loops
  2. write ORCA input  (`n m` header + 0-indexed edge list)
  3. run `orca node 4 in out`  (single-threaded C++)
  4. sum per-node orbit counts (orbits 0-14) and convert to global graphlet
     counts via each orbit's node-multiplicity within its graphlet.

4-node connected graphlet <- orbit (nodes of that orbit per graphlet):
  P4  (4-path)              = sum(orbit4)  / 2
  claw (3-star, K1,3)       = sum(orbit7)  / 1
  C4  (4-cycle)             = sum(orbit8)  / 4
  paw  (triangle+pendant)   = sum(orbit9)  / 1
  diamond (K4 minus edge)   = sum(orbit13) / 2
  K4  (4-clique)            = sum(orbit14) / 4
Also reported (3-node): P3 path = sum(orbit2), triangle = sum(orbit3)/3.

Setup / usage:
  1. build the ORCA binary once:   cd orca_build && make
  2. write per-species ORCA inputs:  python3 orca_count4.py --prep
  3. run ORCA on each input:  for f in orca_build/work/*.in; do
         orca_build/orca node 4 "$f" "${f%.in}.orbits"; done
  4. sum orbits -> results/graphlet4_counts.json:  python3 orca_count4.py --parse
"""

import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ORCA = os.path.join(HERE, "orca_build", "orca")
ROOT = os.path.join(os.path.dirname(HERE), "input")
RESULTS = os.path.join(HERE, "results")


def human(s):
    if s < 90:
        return f"{s:.1f} s"
    if s < 5400:
        return f"{s/60:.1f} min"
    return f"{s/3600:.2f} h"


def convert(tsv_path, orca_in):
    """protein-name edge list -> ORCA integer format; returns (n, m, load_s)."""
    t0 = time.time()
    ids = {}
    edges = []
    seen = set()
    with open(tsv_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == p[1]:
                continue
            a = ids.setdefault(p[0], len(ids))
            b = ids.setdefault(p[1], len(ids))
            key = (a, b) if a < b else (b, a)
            if key not in seen:
                seen.add(key)
                edges.append(key)
    n, m = len(ids), len(edges)
    with open(orca_in, "w") as out:
        out.write(f"{n} {m}\n")
        out.writelines(f"{a} {b}\n" for a, b in edges)
    return n, m, time.time() - t0


def graphlets_from_orbits(orbit_out, n):
    """sum orbit columns 0..14 over all nodes, return global graphlet counts."""
    tot = [0] * 15
    with open(orbit_out) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue
            for j in range(15):
                tot[j] += int(vals[j])
    g = {
        "P4_path":  tot[4] // 2,
        "claw_star": tot[7],
        "C4_cycle": tot[8] // 4,
        "paw":      tot[9],
        "diamond":  tot[13] // 2,
        "K4_clique": tot[14] // 4,
    }
    g["total_connected_4graphlets"] = sum(g.values())
    g["_3node"] = {"P3_path": tot[2], "triangle": tot[3] // 3}
    return g


WORKDIR = os.path.join(HERE, "orca_build", "work")
MANIFEST = os.path.join(WORKDIR, "manifest.json")


def species_dirs(args):
    if args:
        return args
    dirs = sorted(glob.glob(os.path.join(ROOT, "GC*")))
    return [d for d in dirs if os.path.isdir(d)]


def prep(dirs):
    """Trusted step: convert each species tsv -> ORCA .in file; write manifest."""
    os.makedirs(WORKDIR, exist_ok=True)
    manifest = []
    for d in dirs:
        acc = os.path.basename(d.rstrip("/"))
        hits = glob.glob(f"{d.rstrip('/')}/*_network.positive.tsv")
        if not hits:
            print(f"  !! no network file for {acc}")
            continue
        orca_in = os.path.join(WORKDIR, f"{acc}.in")
        orca_out = os.path.join(WORKDIR, f"{acc}.orbits")
        n, m, load_s = convert(hits[0], orca_in)
        print(f"[{acc}] {n:,} nodes  {m:,} edges (avg deg {2*m/n:.1f})  "
              f"-> {os.path.basename(orca_in)} in {human(load_s)}", flush=True)
        manifest.append({"accession": acc, "nodes": n, "edges": m,
                         "in": orca_in, "out": orca_out})
    with open(MANIFEST, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nprepped {len(manifest)} species. ORCA input dir: {WORKDIR}")
    print(f"manifest: {MANIFEST}")


def parse():
    """Trusted step: read .orbits outputs -> graphlet counts -> summary JSON."""
    with open(MANIFEST) as f:
        manifest = json.load(f)
    results = []
    for e in manifest:
        if not os.path.exists(e["out"]):
            print(f"  !! missing ORCA output for {e['accession']} ({e['out']})")
            continue
        g = graphlets_from_orbits(e["out"], e["nodes"])
        results.append({"accession": e["accession"], "nodes": e["nodes"],
                        "edges": e["edges"], "graphlets": g})
        print(f"[{e['accession']}] total connected 4-graphlets = "
              f"{g['total_connected_4graphlets']:,}", flush=True)
    os.makedirs(RESULTS, exist_ok=True)
    out_json = os.path.join(RESULTS, "graphlet4_counts.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_json} ({len(results)} species)")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "--prep":
        prep(species_dirs(sys.argv[2:]))
    elif mode == "--parse":
        parse()
    else:
        sys.exit("usage: orca_count4.py --prep [species_dir ...] | --parse")


if __name__ == "__main__":
    main()
