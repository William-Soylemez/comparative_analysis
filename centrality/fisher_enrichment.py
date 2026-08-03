#!/usr/bin/env python3
"""
GO-term enrichment in the high-betweenness ("connector") protein set, per species.

For each species:
  1. compute exact unweighted betweenness for every network node (cached to
     <acc>_node_betweenness.json so re-runs are instant),
  2. define the high-centrality set H = the top --pct% of nodes by betweenness,
  3. restrict to annotated proteins (universe U = network nodes with >=1 GO term
     in <acc>_GO_map.csv), and
  4. for every GO term appearing in >= --min-count proteins of U, run a one-sided
     Fisher exact test (over-representation: term in H vs term in U\H),
  5. Benjamini-Hochberg FDR-correct across the tested terms.

Design notes:
  * H is the top pct% of ALL network nodes (a purely topological definition),
    then intersected with the annotated universe. So |H| annotated can be a bit
    under pct% of |U| if connector proteins are less annotated -- reported so
    it's transparent. The Fisher test is H vs the rest of the annotated universe.
  * Terms are DIRECT annotations from the GO map by default (consistent with the
    project's earlier go-frequency analyses). --propagate rolls each protein's
    terms up to all is_a/part_of ancestors (excluding the 3 roots) via
    go-basic.obo for a more powerful, more interpretable pass.
  * Term names/namespaces come from go_data/go-basic.obo (reusing goslim_util).

Usage (from this directory):
    ../venv/bin/python fisher_enrichment.py                 # all species, top 5%, direct
    ../venv/bin/python fisher_enrichment.py --pct 5 --propagate
    ../venv/bin/python fisher_enrichment.py --parent .. --min-count 5 --qmax 0.05
"""

import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import fisher_exact
from tqdm import tqdm

from betweenness_threshold_plot import load_graph, norm_betweenness

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from goslim_util import _parse_obo  # noqa: E402

ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}  # BP, MF, CC roots


def node_scores(species_dir, acc, outdir, centrality):
    """Return {node_name: score} for the chosen centrality, using/creating a cache.

    betweenness = exact unweighted normalized betweenness (slow, cached).
    degree      = node degree (cheap; the high set is the most-connected hubs).
    Only the ranking matters for the top-pct% cut, so degree is stored raw.
    """
    cache = os.path.join(outdir, f"{acc}_node_{centrality}.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    g = load_graph(species_dir)
    if centrality == "betweenness":
        vals = norm_betweenness(g)
    elif centrality == "degree":
        vals = g.degree()
    else:
        raise ValueError(f"unknown centrality {centrality!r}")
    d = dict(zip(g.vs["name"], (float(x) for x in vals)))
    with open(cache, "w") as f:
        json.dump(d, f)
    return d


def load_go_map(species_dir, ancestors=None, alt=None, obsolete=None, slim=None):
    """protein -> set(GO terms).

    If `slim` (a GoSlim instance) is given, each protein's direct terms are
    rolled up to their GO-slim categories (union). Else if `ancestors` given,
    terms are propagated to all is_a/part_of ancestors. Otherwise direct terms.
    """
    hits = glob.glob(f"{species_dir.rstrip('/')}/*_GO_map.csv")
    if not hits:
        sys.exit(f"no *_GO_map.csv in {species_dir}")
    import csv
    drop = ROOTS | (obsolete or set())
    m = {}
    with open(hits[0], newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            raw = (row.get("GO_list") or "").strip()
            if not raw:
                continue
            terms = set()
            for t in raw.split(";"):
                t = t.strip()
                if not t:
                    continue
                if alt:
                    t = alt.get(t, t)
                if slim is not None:
                    terms |= slim.map_term(t)      # roll up to slim categories
                else:
                    terms.add(t)
                    if ancestors is not None:
                        terms |= ancestors(t)
            terms -= drop
            if terms:
                m[row["prot_id"]] = terms
    return m


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values for a list of p-values (same order)."""
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1  # BH rank (largest p = rank n)
        val = pvals[idx] * n / i
        prev = min(prev, val)
        q[idx] = prev
    return q


def enrich_species(species_dir, acc, pct, min_count, ancestors, alt, obsolete, slim, names, ns, outdir, centrality):
    score = node_scores(species_dir, acc, outdir, centrality)
    nodes = list(score)
    N = len(nodes)
    k = max(1, math.ceil(pct / 100.0 * N))
    top = set(sorted(nodes, key=lambda n: score[n], reverse=True)[:k])

    go = load_go_map(species_dir, ancestors=ancestors, alt=alt, obsolete=obsolete, slim=slim)
    U = [n for n in nodes if n in go]              # annotated universe
    Uset = set(U)
    H = [n for n in top if n in Uset]              # annotated high-centrality
    nU, nH = len(U), len(H)
    rest = nU - nH

    # term counts in universe and in H
    univ_c, h_c = {}, {}
    for n in U:
        for t in go[n]:
            univ_c[t] = univ_c.get(t, 0) + 1
    for n in H:
        for t in go[n]:
            h_c[t] = h_c.get(t, 0) + 1

    tested = [t for t, c in univ_c.items() if c >= min_count]
    rows, pvals = [], []
    for t in tested:
        a = h_c.get(t, 0)          # in H, has term
        b = nH - a                 # in H, no term
        c = univ_c[t] - a          # in rest, has term
        d = rest - c               # in rest, no term
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        exp = univ_c[t] / nU * nH  # expected count in H under background
        fold = (a / exp) if exp > 0 else float("inf")
        rows.append({
            "term": t, "name": names.get(t, ""), "namespace": ns.get(t, ""),
            "h_count": a, "h_total": nH, "univ_count": univ_c[t], "univ_total": nU,
            "fold_enrichment": round(fold, 2),
            "odds_ratio": (round(odds, 3) if math.isfinite(odds) else None),
            "p": p,
        })
        pvals.append(p)

    if rows:
        q = bh_fdr(pvals)
        for row, qi in zip(rows, q):
            row["q"] = float(qi)
        rows.sort(key=lambda r: (r["q"], r["p"]))
    return {
        "species": acc, "n_nodes": N, "pct": pct, "k_top_all": k,
        "n_universe": nU, "n_high_annotated": nH, "n_terms_tested": len(tested),
        "rows": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parent", default="../input")
    ap.add_argument("--centrality", choices=["betweenness", "degree"], default="betweenness",
                    help="node ranking for the high-centrality set (default betweenness)")
    ap.add_argument("--pct", type=float, default=5.0, help="top %% by centrality (default 5)")
    ap.add_argument("--min-count", type=int, default=5,
                    help="min proteins in universe annotated with a term to test it (default 5)")
    ap.add_argument("--qmax", type=float, default=0.05, help="FDR threshold for the 'significant' list")
    ap.add_argument("--propagate", action="store_true", help="roll terms up to GO ancestors")
    ap.add_argument("--slim", nargs="?", const="../go_data/goslim_generic.obo", default=None,
                    help="roll terms up to GO-slim categories (default slim: goslim_generic.obo)")
    ap.add_argument("--obo", default="../go_data/go-basic.obo")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    species_dirs = [d for d in sorted(glob.glob(os.path.join(args.parent, "GC[AF]_*")))
                    if os.path.isdir(d) and glob.glob(os.path.join(d, "*_network.positive.tsv"))]
    if not species_dirs:
        raise SystemExit(f"no species dirs under {args.parent}")

    print("parsing go-basic.obo ...", flush=True)
    t0 = time.time()
    names, ns, parents, obsolete, alt = _parse_obo(args.obo)
    print(f"  {len(names):,} terms ({time.time()-t0:.1f}s)")

    anc_cache = {}
    def ancestors(term):
        if term in anc_cache:
            return anc_cache[term]
        seen, stack = set(), [term]
        while stack:
            x = stack.pop()
            for p in parents.get(x, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        anc_cache[term] = seen
        return seen
    anc = ancestors if args.propagate else None

    slim = None
    if args.slim:
        from goslim_util import GoSlim
        slim = GoSlim(args.obo, args.slim)
        print(f"  GO-slim: {len(slim.slim)} categories from {os.path.basename(args.slim)}")

    base = "slim" if args.slim else ("propagated" if args.propagate else "direct")
    # betweenness stays unprefixed (backward compat); degree is prefixed
    tag = base if args.centrality == "betweenness" else f"{args.centrality}_{base}"
    summary = {}
    for sdir in tqdm(species_dirs, desc="enrichment", unit="sp"):
        acc = os.path.basename(sdir.rstrip("/"))
        res = enrich_species(sdir, acc, args.pct, args.min_count, anc, alt, obsolete, slim, names, ns, args.outdir, args.centrality)
        sig = [r for r in res["rows"] if r.get("q", 1.0) <= args.qmax]

        # per-species TSV of all tested terms (sorted by q)
        tsv = os.path.join(args.outdir, f"{acc}_enrichment_{tag}.tsv")
        with open(tsv, "w") as f:
            f.write("term\tname\tnamespace\th_count\th_total\tuniv_count\tuniv_total"
                    "\tfold_enrichment\todds_ratio\tp\tq\n")
            for r in res["rows"]:
                f.write(f"{r['term']}\t{r['name']}\t{r['namespace']}\t{r['h_count']}"
                        f"\t{r['h_total']}\t{r['univ_count']}\t{r['univ_total']}"
                        f"\t{r['fold_enrichment']}\t{r['odds_ratio']}\t{r['p']:.3e}"
                        f"\t{r['q']:.3e}\n")
        summary[acc] = {
            "n_nodes": res["n_nodes"], "n_universe": res["n_universe"],
            "n_high_annotated": res["n_high_annotated"],
            "n_terms_tested": res["n_terms_tested"], "n_significant": len(sig),
            "top_significant": [
                {"term": r["term"], "name": r["name"], "namespace": r["namespace"],
                 "h_count": r["h_count"], "univ_count": r["univ_count"],
                 "fold_enrichment": r["fold_enrichment"], "q": round(r["q"], 4)}
                for r in sig[:15]
            ],
        }
        tqdm.write(f"  {acc}: U={res['n_universe']} H={res['n_high_annotated']} "
                   f"tested={res['n_terms_tested']}  significant(q<={args.qmax})={len(sig)}")

    outjson = os.path.join(args.outdir, f"enrichment_summary_{tag}.json")
    with open(outjson, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote {outjson}")
    print(f"per-species: <acc>_enrichment_{tag}.tsv  (+ <acc>_node_{args.centrality}.json cache)")


if __name__ == "__main__":
    main()
