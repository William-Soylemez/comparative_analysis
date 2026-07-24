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
"""

import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ORCA = os.path.join(HERE, "orca_build", "orca")
ROOT = os.path.dirname(HERE)


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


def run_species(species_dir, workdir):
    acc = os.path.basename(species_dir.rstrip("/"))
    hits = glob.glob(f"{species_dir.rstrip('/')}/*_network.positive.tsv")
    if not hits:
        print(f"  !! no network file for {acc}")
        return None
    orca_in = os.path.join(workdir, f"{acc}.in")
    orca_out = os.path.join(workdir, f"{acc}.orbits")
    n, m, load_s = convert(hits[0], orca_in)
    print(f"[{acc}] {n:,} nodes  {m:,} edges (avg deg {2*m/n:.1f})  "
          f"converted in {human(load_s)}", flush=True)
    t = time.time()
    r = subprocess.run([ORCA, "node", "4", orca_in, orca_out],
                       capture_output=True, text=True)
    run_s = time.time() - t
    if r.returncode != 0 or not os.path.exists(orca_out):
        print(f"  !! ORCA failed: {r.stderr.strip() or r.stdout.strip()}")
        return None
    g = graphlets_from_orbits(orca_out, n)
    print(f"[{acc}] ORCA done in {human(run_s)}  ->  "
          f"total connected 4-graphlets = {g['total_connected_4graphlets']:,}",
          flush=True)
    return {"accession": acc, "nodes": n, "edges": m,
            "orca_seconds": round(run_s, 2), "graphlets": g}


def main():
    workdir = os.path.join(HERE, "orca_build", "work")
    os.makedirs(workdir, exist_ok=True)
    if len(sys.argv) > 1:
        dirs = sys.argv[1:]
    else:
        dirs = sorted(glob.glob(os.path.join(ROOT, "GC*")))
        dirs = [d for d in dirs if os.path.isdir(d)]
    results = []
    t0 = time.time()
    for d in dirs:
        res = run_species(d, workdir)
        if res:
            results.append(res)
    out_json = os.path.join(HERE, "graphlet4_counts.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nAll done in {human(time.time()-t0)}. Wrote {out_json}")


if __name__ == "__main__":
    main()
