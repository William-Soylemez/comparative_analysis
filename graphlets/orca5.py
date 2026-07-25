#!/usr/bin/env python3
"""
n=5 connected-graphlet counting via ORCA `node 5` (73 per-node orbits, 0-72).
Orbits 15-72 belong to the 21 connected 5-node graphlets.

We do NOT hardcode the orbit->graphlet reduction. Instead we self-calibrate:
run ORCA on each isolated 5-node graphlet H (where count(H)=1 by construction),
so each orbit column-sum over H's 5 nodes equals that orbit's node-multiplicity
within its graphlet. Then on any real graph, count(Gk) = colsum(rep_orbit_k)/mult.

Modes:
  --gen        generate the 21 connected 5-node graphlets + ORCA inputs (calib/)
  --table      read calib ORCA outputs -> reduction table (orbit5_reduction.json)
  --validate   brute-force connected-5-subgraph counts on a small random graph,
               compare to ORCA-derived counts (needs the table + orca binary run
               on the validation graph, done via bash)
  --parse a b  parse real-species .orbits5 files for given accessions -> counts
"""
import glob
import json
import os
import sys
from itertools import combinations

import igraph as ig

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "orca_build")
CALIB = os.path.join(BUILD, "calib")
WORK = os.path.join(BUILD, "work")
TABLE = os.path.join(HERE, "orbit5_reduction.json")
ALL_PAIRS = list(combinations(range(5), 2))


def canon_cert(edges, n=5):
    g = ig.Graph(n=n, edges=edges, directed=False)
    perm = g.canonical_permutation()
    g = g.permute_vertices(perm)
    return tuple(sorted((min(a, b), max(a, b)) for a, b in
                        (e.tuple for e in g.es)))


def gen():
    """enumerate the 21 connected labelled-up-to-iso 5-node graphs."""
    os.makedirs(CALIB, exist_ok=True)
    reps = {}  # cert -> representative edge list
    for mask in range(1 << len(ALL_PAIRS)):
        edges = [ALL_PAIRS[i] for i in range(len(ALL_PAIRS)) if (mask >> i) & 1]
        if len(edges) < 4:  # need >=4 edges to be connected on 5 nodes
            continue
        g = ig.Graph(n=5, edges=edges, directed=False)
        if not g.is_connected():
            continue
        cert = canon_cert(edges)
        if cert not in reps:
            reps[cert] = edges
    # deterministic order: by (#edges, sorted degree sequence, cert)
    items = []
    for cert, edges in reps.items():
        g = ig.Graph(n=5, edges=edges, directed=False)
        items.append((len(edges), tuple(sorted(g.degree(), reverse=True)),
                      cert, edges))
    items.sort(key=lambda t: (t[0], t[1], t[2]))
    meta = []
    for gid, (ne, degseq, cert, edges) in enumerate(items):
        fn = os.path.join(CALIB, f"gl{gid:02d}.in")
        with open(fn, "w") as f:
            f.write(f"5 {len(edges)}\n")
            f.writelines(f"{a} {b}\n" for a, b in edges)
        meta.append({"gid": gid, "edges": ne, "degseq": list(degseq),
                     "cert": [list(e) for e in cert], "name": name_graphlet(ne, degseq)})
    with open(os.path.join(CALIB, "graphlets.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"generated {len(meta)} connected 5-node graphlets in {CALIB}")
    for m in meta:
        print(f"  gl{m['gid']:02d}: {m['edges']} edges, deg {m['degseq']}  {m['name']}")


def name_graphlet(ne, degseq):
    d = tuple(sorted(degseq, reverse=True))
    known = {
        (4, (2, 2, 2, 1, 1)): "P5 (5-path)",
        (4, (4, 1, 1, 1, 1)): "S4 (4-star/claw)",
        (4, (3, 2, 1, 1, 1)): "chair/fork",
        (5, (2, 2, 2, 2, 2)): "C5 (5-cycle)",
        (10, (4, 4, 4, 4, 4)): "K5 (5-clique)",
    }
    return known.get((ne, d), f"{ne}e_{'.'.join(map(str, d))}")


def read_orbit_sums(path, ncols=73):
    tot = [0] * ncols
    with open(path) as f:
        for line in f:
            v = line.split()
            if not v:
                continue
            for j in range(ncols):
                tot[j] += int(v[j])
    return tot


def build_table():
    meta = json.load(open(os.path.join(CALIB, "graphlets.json")))
    reduction = {}   # gid -> {"rep_orbit": o, "mult": mult, "orbits": {o: mult}}
    orbit_owner = {}  # orbit(15-72) -> gid  (sanity: each owned by one graphlet)
    for m in meta:
        gid = m["gid"]
        sums = read_orbit_sums(os.path.join(CALIB, f"gl{gid:02d}.orbits5"))
        # in isolated H, count(H)=1 so colsum(o)=mult(o) for o belonging to H
        orbits = {o: sums[o] for o in range(15, 73) if sums[o] > 0}
        if not orbits:
            raise SystemExit(f"gl{gid:02d}: no 5-node orbits fired?!")
        # consistency: multiplicities must sum to 5 (partition of the 5 nodes)
        if sum(orbits.values()) != 5:
            raise SystemExit(f"gl{gid:02d}: orbit mults {orbits} don't sum to 5")
        for o in orbits:
            if o in orbit_owner:
                raise SystemExit(f"orbit {o} claimed by gl{orbit_owner[o]} and gl{gid}")
            orbit_owner[o] = gid
        rep = min(orbits)
        reduction[str(gid)] = {"rep_orbit": rep, "mult": orbits[rep],
                               "orbits": {str(o): mm for o, mm in orbits.items()},
                               "edges": m["edges"], "degseq": m["degseq"],
                               "name": m["name"]}
    json.dump(reduction, open(TABLE, "w"), indent=2)
    print(f"built reduction table for {len(reduction)} graphlets -> {TABLE}")
    print(f"orbits 15-72 covered: {len(orbit_owner)} (expect 58)")


def counts_from_orbits(path):
    """real graph .orbits5 -> {gid: count}, using self-consistency across orbits."""
    red = json.load(open(TABLE))
    sums = read_orbit_sums(path)
    out = {}
    for gid, r in red.items():
        vals = set()
        for o, mult in r["orbits"].items():
            vals.add(sums[int(o)] // mult)
            if sums[int(o)] % mult != 0:
                vals.add(("noninteger", o, sums[int(o)], mult))
        cnt = sums[r["rep_orbit"]] // r["mult"]
        # every orbit of the graphlet must give the same count (built-in check)
        clean = {v for v in vals if isinstance(v, int)}
        out[gid] = {"count": cnt, "name": r["name"], "edges": r["edges"],
                    "consistent": len(clean) == 1 and not any(
                        not isinstance(v, int) for v in vals)}
    return out


def parse(accessions):
    red = json.load(open(TABLE))
    results = []
    for acc in accessions:
        path = os.path.join(WORK, f"{acc}.orbits5")
        if not os.path.exists(path):
            print(f"  !! missing {path}")
            continue
        c = counts_from_orbits(path)
        total = sum(v["count"] for v in c.values())
        allok = all(v["consistent"] for v in c.values())
        results.append({"accession": acc, "total_connected_5graphlets": total,
                        "all_orbits_consistent": allok,
                        "by_class": {gid: c[gid] for gid in sorted(c, key=int)}})
        print(f"[{acc}] total 5-graphlets = {total:,}  (orbit-consistency: {allok})")
    out = os.path.join(HERE, "graphlet5_counts.json")
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "--gen":
        gen()
    elif mode == "--table":
        build_table()
    elif mode == "--parse":
        parse(sys.argv[2:])
    else:
        sys.exit("usage: orca5.py --gen | --table | --parse <acc...>")
