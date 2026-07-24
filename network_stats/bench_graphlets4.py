#!/usr/bin/env python3
"""
Benchmark n=4 connected-graphlet (motif) counting on one species' PPI network,
single-threaded, using igraph's C RAND-ESU implementation.

Two modes:
  --project : load graph, run ESU over a SAMPLE of root vertices to project the
              full exhaustive single-core runtime and estimate the total count.
              Cheap; used to decide whether the real run is worth kicking off.
  --full    : run the exhaustive exact size-4 motif count and report wall time.

igraph.motifs_randesu(size=4) returns per-isomorphism-class counts of connected
INDUCED 4-node subgraphs (the 6 connected 4-graphlets; disconnected classes are
reported as NaN). motifs_randesu_estimate does the same over a vertex sample.
"""

import argparse
import glob
import os
import sys
import time

import igraph as ig


def human(seconds):
    if seconds < 90:
        return f"{seconds:.1f} s"
    if seconds < 5400:
        return f"{seconds/60:.1f} min"
    if seconds < 172800:
        return f"{seconds/3600:.2f} h"
    return f"{seconds/86400:.1f} days"


def load_graph(species_dir):
    hits = glob.glob(f"{species_dir.rstrip('/')}/*_network.positive.tsv")
    if not hits:
        sys.exit(f"no *_network.positive.tsv in {species_dir}")
    t0 = time.time()
    edges = []
    seen = set()
    verts = {}
    def vid(name):
        i = verts.get(name)
        if i is None:
            i = len(verts)
            verts[name] = i
        return i
    with open(hits[0]) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != p[1]:
                a, b = vid(p[0]), vid(p[1])
                key = (a, b) if a < b else (b, a)
                if key not in seen:
                    seen.add(key)
                    edges.append(key)
    g = ig.Graph(n=len(verts), edges=edges, directed=False)
    g.simplify()  # belt-and-suspenders: drop any residual multi/self edges
    load_t = time.time() - t0
    return g, hits[0], load_t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species_dir")
    ap.add_argument("--project", action="store_true", help="sampled projection only")
    ap.add_argument("--full", action="store_true", help="exhaustive full count")
    ap.add_argument("--sample", type=int, default=200,
                    help="root vertices to sample in projection mode (default 200)")
    ap.add_argument("--size", type=int, default=4)
    args = ap.parse_args()

    g, path, load_t = load_graph(args.species_dir)
    n, m = g.vcount(), g.ecount()
    acc = os.path.basename(args.species_dir.rstrip("/"))
    print(f"species:   {acc}", flush=True)
    print(f"file:      {os.path.basename(path)}", flush=True)
    print(f"graph:     {n:,} nodes, {m:,} edges  (avg degree {2*m/n:.1f})", flush=True)
    print(f"load time: {human(load_t)}", flush=True)
    print(f"connected components: {len(g.connected_components())}", flush=True)

    # size-3 is cheap; useful sanity signal and a lower bound on cost
    t = time.time()
    m3 = g.motifs_randesu(size=3)
    t3 = time.time() - t
    tri = m3[3] if len(m3) > 3 and m3[3] == m3[3] else float("nan")
    print(f"\nsize-3 motifs (exact): {t3:.2f} s  -> classes {m3}", flush=True)

    if args.project:
        frac = args.sample / n
        print(f"\n--- projection: sampling {args.sample} of {n:,} root vertices "
              f"({100*frac:.2f}%) ---", flush=True)
        t = time.time()
        est = g.motifs_randesu_estimate(size=args.size, sample=args.sample)
        samp_t = time.time() - t
        # exhaustive time ~ sample_time / fraction_of_roots_visited
        proj = samp_t / frac
        print(f"sampled ESU wall:      {human(samp_t)} for {args.sample} roots", flush=True)
        print(f"estimated total count: {est:,.0f}  (size-{args.size} connected graphlets)", flush=True)
        print(f"\nPROJECTED exhaustive single-core runtime: ~{human(proj)}", flush=True)
        print(f"  (linear extrapolation from {100*frac:.2f}% of root vertices; "
              f"order-of-magnitude only)", flush=True)

    if args.full:
        print(f"\n--- exhaustive exact size-{args.size} count (single-threaded) ---",
              flush=True)
        t = time.time()
        counts = g.motifs_randesu(size=args.size)
        full_t = time.time() - t
        total = sum(c for c in counts if c == c)  # skip NaN (disconnected classes)
        print(f"WALL TIME: {human(full_t)}  ({full_t:.1f} s)", flush=True)
        print(f"total connected size-{args.size} graphlets: {total:,.0f}", flush=True)
        print(f"per-isomorphism-class counts: {counts}", flush=True)


if __name__ == "__main__":
    main()
