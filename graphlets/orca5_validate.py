#!/usr/bin/env python3
"""
Independent validation of the n=5 ORCA reduction: brute-force count connected
INDUCED 5-node subgraphs on a small random graph, bin by canonical form, and
compare per-class to the ORCA-derived counts.

  --gen    : build a small random graph, write ORCA input (valid.in) and the
             brute-force per-graphlet counts (valid_brute.json)
  --check  : parse valid.orbits5 (produced by `orca node 5` between the two)
             with the reduction table, compare to brute force.
"""
import json
import os
import sys
from itertools import combinations

import igraph as ig

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "orca_build")
CALIB = os.path.join(BUILD, "calib")
VIN = os.path.join(BUILD, "valid.in")
VORB = os.path.join(BUILD, "valid.orbits5")
VBRUTE = os.path.join(BUILD, "valid_brute.json")
TABLE = os.path.join(HERE, "results", "orbit5_reduction.json")


def cert_of(g):
    perm = g.canonical_permutation()
    g = g.permute_vertices(perm)
    return tuple(sorted((min(a, b), max(a, b)) for a, b in
                        (e.tuple for e in g.es)))


def gen(n=30, p=0.30, seed=7):
    import random
    random.seed(seed)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)
             if random.random() < p]
    G = ig.Graph(n=n, edges=edges, directed=False)
    with open(VIN, "w") as f:
        f.write(f"{n} {len(edges)}\n")
        f.writelines(f"{a} {b}\n" for a, b in edges)
    # cert -> gid map from calibration graphlets
    meta = json.load(open(os.path.join(CALIB, "graphlets.json")))
    cert2gid = {tuple(tuple(e) for e in m["cert"]): m["gid"] for m in meta}
    # brute force over all 5-subsets
    counts = {m["gid"]: 0 for m in meta}
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    checked = 0
    for combo in combinations(range(n), 5):
        s = set(combo)
        # induced edges
        ie = [(x, y) for x, y in combinations(combo, 2) if y in adj[x]]
        if len(ie) < 4:
            continue
        idx = {v: k for k, v in enumerate(combo)}
        h = ig.Graph(n=5, edges=[(idx[x], idx[y]) for x, y in ie], directed=False)
        if not h.is_connected():
            continue
        cert = cert_of(h)
        gid = cert2gid.get(cert)
        if gid is None:
            raise SystemExit(f"unknown cert for subset {combo}: {cert}")
        counts[gid] += 1
        checked += 1
    json.dump({"n": n, "edges": len(edges), "counts": counts,
               "total": checked}, open(VBRUTE, "w"), indent=2)
    print(f"validation graph: {n} nodes, {len(edges)} edges")
    print(f"brute-force connected induced 5-subgraphs: {checked:,}")


def check():
    brute = json.load(open(VBRUTE))
    red = json.load(open(TABLE))
    # ORCA-derived from valid.orbits5
    tot = [0] * 73
    for line in open(VORB):
        v = line.split()
        if v:
            for j in range(73):
                tot[j] += int(v[j])
    ok = True
    print(f"{'gid':>3} {'name':28} {'brute':>12} {'ORCA':>12}  match")
    for gid_s, r in sorted(red.items(), key=lambda kv: int(kv[0])):
        gid = int(gid_s)
        orca = tot[r["rep_orbit"]] // r["mult"]
        b = brute["counts"][str(gid)] if str(gid) in brute["counts"] else brute["counts"][gid]
        m = (orca == b)
        ok = ok and m
        print(f"{gid:>3} {r['name'][:28]:28} {b:>12,} {orca:>12,}  {'OK' if m else 'XXXX'}")
    tb = brute["total"]
    to = sum((tot[r["rep_orbit"]] // r["mult"]) for r in red.values())
    print(f"\ntotal brute={tb:,}  ORCA={to:,}  {'ALL MATCH' if ok and tb==to else 'MISMATCH'}")


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else ""
    if m == "--gen":
        gen()
    elif m == "--check":
        check()
    else:
        sys.exit("usage: orca5_validate.py --gen | --check")
