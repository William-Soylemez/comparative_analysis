#!/usr/bin/env python3
"""
Benchmark whether exact/average node connectivity is realistic for one species,
and whether multiprocessing + pair-sampling makes it viable.

average_node_connectivity(G) = mean over all C(n,2) node pairs of the local
node connectivity kappa(u,v) (a max-flow per pair). It is:
  * embarrassingly parallel ACROSS pairs (each kappa is independent), but
  * pure-Python in networkx and CPU-bound, so real speedup needs PROCESSES
    (multiprocessing), not threads -- the GIL blocks thread parallelism here.

This script measures single-core per-pair cost, measures parallel throughput
across --procs processes, projects the full exact runtime, and then shows the
sampling estimate (mean kappa over K random pairs +/- standard error), which is
what actually makes this tractable.

Usage (on a compute node, NOT the login node):
    python3 bench_node_connectivity.py SPECIES_DIR [--pairs 500] [--procs N] [--serial-sample 15]

Example:
    python3 bench_node_connectivity.py ../input/GCF_000182965.3 --pairs 500
    # or grab all cores on the node:
    python3 bench_node_connectivity.py ../input/GCF_000182965.3 --pairs 2000 --procs 144
"""

import argparse
import glob
import math
import os
import random
import sys
import time
from statistics import mean, pstdev

import networkx as nx
from networkx.algorithms.connectivity import build_auxiliary_node_connectivity
from networkx.algorithms.connectivity.connectivity import local_node_connectivity
from networkx.algorithms.flow import build_residual_network

# module globals so forked workers inherit the graph without pickling it
_G = None
_H = None
_R = None


def load_largest_cc(species_dir):
    hits = glob.glob(f"{species_dir.rstrip('/')}/*_network.positive.tsv")
    if not hits:
        sys.exit(f"no *_network.positive.tsv in {species_dir}")
    g = nx.Graph()
    with open(hits[0]) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] != p[1]:
                g.add_edge(p[0], p[1])
    cc = max(nx.connected_components(g), key=len)
    return g.subgraph(cc).copy()


def _init():
    # one auxiliary digraph + residual network per worker, reused across pairs
    global _H, _R
    _H = build_auxiliary_node_connectivity(_G)
    _R = build_residual_network(_H, "capacity")


def _kappa(pair):
    u, v = pair
    return local_node_connectivity(_G, u, v, auxiliary=_H, residual=_R)


def human(seconds):
    if seconds < 90:
        return f"{seconds:.1f} s"
    if seconds < 5400:
        return f"{seconds/60:.1f} min"
    if seconds < 172800:
        return f"{seconds/3600:.1f} h"
    return f"{seconds/86400:.1f} days"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species_dir")
    ap.add_argument("--pairs", type=int, default=500, help="random pairs to sample for the estimate (default 500)")
    ap.add_argument("--procs", type=int, default=os.cpu_count(), help="worker processes (default: all cores)")
    ap.add_argument("--serial-sample", type=int, default=15, help="pairs to time on a single core for per-pair cost")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    random.seed(args.seed)

    global _G
    t0 = time.time()
    _G = load_largest_cc(args.species_dir)
    n, m = _G.number_of_nodes(), _G.number_of_edges()
    total_pairs = n * (n - 1) // 2
    acc = os.path.basename(args.species_dir.rstrip("/"))
    print(f"species: {acc}")
    print(f"largest CC: {n:,} nodes, {m:,} edges  (loaded in {human(time.time()-t0)})")
    print(f"total pairs C(n,2): {total_pairs:,}")
    nodes = list(_G)

    pairs = [tuple(random.sample(nodes, 2)) for _ in range(args.pairs)]

    # ---- single-core per-pair cost (time-guarded) ----
    _init()  # build aux/residual in the main process
    ns = min(args.serial_sample, len(pairs))
    t = time.time()
    done = 0
    for pr in pairs[:ns]:
        _kappa(pr)
        done += 1
        if time.time() - t > 30:  # bail if a single core is very slow
            break
    per_pair = (time.time() - t) / done
    print(f"\n--- single core ---")
    print(f"timed {done} pairs: {per_pair*1000:.1f} ms/pair")
    print(f"exact full run, 1 core: {total_pairs:,} pairs -> {human(total_pairs*per_pair)}")

    # ---- parallel throughput across --procs processes ----
    print(f"\n--- {args.procs} processes ---")
    t = time.time()
    # fork so workers inherit _G copy-on-write (no pickling the graph 144x);
    # TACC is Linux where fork is the default anyway.
    from multiprocessing import get_context
    ctx = get_context("fork")
    with ctx.Pool(processes=args.procs, initializer=_init) as pool:
        kappas = pool.map(_kappa, pairs, chunksize=max(1, len(pairs) // (args.procs * 4)))
    wall = time.time() - t
    thru = len(pairs) / wall
    print(f"{len(pairs):,} pairs in {human(wall)}  ({thru:,.0f} pairs/s, {args.procs*wall/len(pairs)*1000:.1f} ms/pair/core incl. overhead)")
    print(f"exact full run, {args.procs} cores: {human(total_pairs/thru)}")

    # ---- sampling estimate (the practical option) ----
    mu, sd = mean(kappas), pstdev(kappas)
    se = sd / math.sqrt(len(kappas))
    ci = 1.96 * se
    print(f"\n--- sampling estimate ---")
    print(f"avg node connectivity ~ {mu:.3f}  +/- {ci:.3f}  (95% CI, from {len(kappas):,} pairs)")
    print(f"relative 95% half-width: {100*ci/mu:.2f}%" if mu else "mean is 0")
    if mu and sd:
        need = math.ceil((1.96 * sd / (0.02 * mu)) ** 2)  # pairs for +/-2%
        print(f"pairs needed for +/-2%: ~{need:,}  ->  {human(need/thru)} on {args.procs} cores")
    print(f"\nthis {len(pairs):,}-pair estimate cost {human(wall)} on {args.procs} cores.")


if __name__ == "__main__":
    main()
