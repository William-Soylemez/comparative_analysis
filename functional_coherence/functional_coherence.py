#!/usr/bin/env python3
"""
Functional coherence — PHILHARMONIC-style per-cluster GO agreement.

Reproduces the cluster-coherence metric from PHILHARMONIC's
nb/02_functional_permutation_analysis.ipynb: for each cluster, take all pairs of
member proteins (the cluster's `members` PLUS the ReCIPE re-added proteins at the
0.75 degree threshold) and compute the Jaccard similarity of their GO-term sets
(|A n B| / |A u B|). The *mean* over pairs is the cluster's coherence score.
Collecting one score per cluster gives a per-species distribution; overlaying
those across the fungal species shows how functionally clean each species'
modules are.

NB on the metric: the reference uses `1 - scipy.spatial.distance.jaccard(bv1,
bv2)` on GO bit vectors, and MEAN (not median) over pairs. Two facts about it:
  - Median is useless here: >80% of within-cluster pairs share zero GO terms, so
    the per-cluster median is 0 for ~88% of clusters. Mean is the informative
    (and canonical) summary.
  - A pair of two UNANNOTATED proteins scores 1.0 (scipy treats two all-zero bit
    vectors as identical). So low-annotation species get coherence spuriously
    *inflated*. This is the annotation-depth confound, and it runs opposite to
    intuition.

We therefore compute two variants per cluster:
  phil       : faithful to the reference (mean; both-unannotated pair = 1.0).
  annotated  : proteins with no GO terms are dropped before pairing, so the
               "two blanks look identical" inflation can't happen. Cleaner for
               cross-species comparison, at the cost of departing from the ref.

Per-species GO-annotation coverage (fraction of clustered proteins with >=1 GO
term) is reported alongside so the confound is visible.

Usage:
    python3 functional_coherence.py [SPECIES_DIR ...]

With no args, discovers every ../<acc>/ dir containing a *_clusters.json.
Writes functional_coherence.json and functional_coherence.png here.
"""

import sys
import os
import csv
import json
import glob
from itertools import combinations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RECIPE_THRESHOLD = "0.75"  # ReCIPE degree threshold used by the reference notebook
MIN_CLUSTER = 2            # need >= this many proteins to form >=1 pair

# Accession -> (short organism label, phylum) for nicer plots. Order here is the
# rough phylogeny (early-diverging -> crown) used to sort the plot.
SPECIES_META = {
    "GCF_000203795.2": ("B. dendrobatidis", "Chytridiomycota"),
    "GCA_025594325.1": ("B. emersonii", "Blastocladiomycota"),
    "GCF_025094135.1": ("M. mucedo", "Mucoromycota"),
    "GCF_026210795.1": ("R. irregularis", "Mucoromycota"),
    "GCF_025024165.1": ("K. alabastrina", "Zoopagomycota"),
    "GCF_000146465.1": ("E. intestinalis", "Microsporidia"),
    "GCF_000146045.2": ("S. cerevisiae", "Ascomycota"),
    "GCF_000182965.3": ("C. albicans", "Ascomycota"),
    "GCF_000182925.2": ("N. crassa", "Ascomycota"),
    "GCF_000149245.1": ("C. neoformans", "Basidiomycota"),
    "GCA_014872705.1": ("A. bisporus", "Basidiomycota"),
}


def discover_species(args):
    if args:
        dirs = [a.rstrip("/") for a in args]
    else:
        # look one level up (script lives in functional_coherence/)
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        dirs = sorted(
            os.path.join(root, d) for d in os.listdir(root)
            if os.path.isdir(os.path.join(root, d))
            and glob.glob(f"{os.path.join(root, d)}/*_clusters.json")
        )
    out = []
    for d in dirs:
        hits = glob.glob(f"{d}/*_clusters.json")
        if not hits:
            print(f"  skip {d}: no *_clusters.json", file=sys.stderr)
            continue
        acc = os.path.basename(hits[0]).replace("_clusters.json", "")
        out.append((acc, d))
    return out


def load_go_map(go_map_path):
    """prot_id -> frozenset of GO ids (empty frozenset if unannotated)."""
    p2go = {}
    with open(go_map_path) as f:
        for row in csv.DictReader(f):
            terms = frozenset(g for g in (row.get("GO_list") or "").split(";") if g)
            p2go[row["prot_id"]] = terms
    return p2go


def pair_jaccard(a, b):
    """1 - scipy.jaccard on GO bit vectors, computed set-wise (identical result).
    Two unannotated proteins (empty|empty) -> 1.0, matching the reference."""
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def cluster_members(clust):
    """The protein set the reference pairs over: members + ReCIPE 0.75 re-adds."""
    recipe = clust.get("recipe", {}).get("degree", {}).get(RECIPE_THRESHOLD, [])
    # de-dup while preserving membership; recipe adds are disjoint from members in practice
    return list(dict.fromkeys(list(clust["members"]) + list(recipe)))


def cluster_coherence(members, p2go, mode):
    """Mean pairwise Jaccard over a cluster. Returns score or None if < 1 pair.

    mode == 'phil'      : faithful to reference (both-unannotated pair = 1.0).
    mode == 'annotated' : drop members with empty GO set before pairing.
    """
    sets = [p2go.get(m, frozenset()) for m in members]
    if mode == "annotated":
        sets = [s for s in sets if s]
    if len(sets) < 2:
        return None
    jac = [pair_jaccard(a, b) for a, b in combinations(sets, 2)]
    return float(np.mean(jac)) if jac else None


def score_species(acc, d):
    clusters = json.load(open(f"{d}/{acc}_clusters.json"))
    p2go = load_go_map(f"{d}/{acc}_GO_map.csv")

    per_cluster = {"phil": [], "annotated": []}
    clustered = set()
    for c in clusters.values():
        members = cluster_members(c)
        clustered.update(members)
        if len(members) < MIN_CLUSTER:
            continue
        for mode in ("phil", "annotated"):
            score = cluster_coherence(members, p2go, mode)
            if score is not None:
                per_cluster[mode].append(score)

    annotated = sum(1 for m in clustered if p2go.get(m))
    coverage = annotated / len(clustered) if clustered else 0.0

    return {
        "n_clusters_total": len(clusters),
        "n_scored_phil": len(per_cluster["phil"]),
        "n_clustered_proteins": len(clustered),
        "annotation_coverage": coverage,
        "scores_phil": per_cluster["phil"],
        "scores_annotated": per_cluster["annotated"],
    }


def _panel(ax, data, cov, labels, colors, title, ylabel, show_cov):
    parts = ax.violinplot(data, showextrema=False, widths=0.85)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[i])
        body.set_edgecolor(colors[i])
        body.set_alpha(0.55)
    ax.boxplot(data, widths=0.13, showfliers=False, patch_artist=True,
               medianprops=dict(color="#222222", linewidth=1.3),
               whiskerprops=dict(color="#777777"),
               capprops=dict(color="#777777"),
               boxprops=dict(facecolor="white", edgecolor="#777777"))
    ymax = max((max(d) for d in data if d), default=1.0)
    for i, vals in enumerate(data):
        ax.text(i + 1, ymax * 1.04, f"{np.mean(vals):.3f}", ha="center",
                fontsize=7.5, color="#222222")
        if show_cov:
            ax.text(i + 1, -ymax * 0.11, f"cov {cov[i]:.0%}", ha="center",
                    fontsize=7, color="#b06a00")
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=9, fontstyle="italic")
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_ylim(-ymax * 0.16, ymax * 1.12)
    ax.grid(True, axis="y", color="#e6e6e3", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=11)


def plot(results, order, outfile):
    labels = [SPECIES_META.get(a, (a, ""))[0] for a in order]
    cov = [results[a]["annotation_coverage"] for a in order]
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(1, len(order) - 1)) for i in range(len(order))]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 11))
    _panel(ax1, [results[a]["scores_phil"] for a in order], cov, labels, colors,
           "Faithful PHILHARMONIC coherence — mean pairwise GO Jaccard per cluster\n"
           "(numbers = mean over clusters; two-unannotated pairs count as 1.0, so low-coverage "
           "species are inflated — see coverage)",
           "Cluster coherence (phil)", show_cov=True)
    _panel(ax2, [results[a]["scores_annotated"] for a in order], cov, labels, colors,
           "Annotated-only variant — unannotated proteins dropped before pairing\n"
           "(removes the two-blanks-look-identical inflation; cleaner cross-species comparison)",
           "Cluster coherence (annotated only)", show_cov=False)
    fig.suptitle("Functional coherence of PHILHARMONIC clusters across 11 fungi",
                 fontsize=13, y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  wrote {outfile}")


def main():
    species = discover_species(sys.argv[1:])
    if not species:
        print("No *_clusters.json found.", file=sys.stderr)
        sys.exit(1)

    results = {}
    for acc, d in species:
        r = score_species(acc, d)
        results[acc] = r
        m_phil = np.mean(r["scores_phil"]) if r["scores_phil"] else float("nan")
        m_ann = np.mean(r["scores_annotated"]) if r["scores_annotated"] else float("nan")
        print(f"  {SPECIES_META.get(acc, (acc,''))[0]:<18} "
              f"clusters={r['n_scored_phil']:>4}/{r['n_clusters_total']:<4} "
              f"coverage={r['annotation_coverage']:.0%}  "
              f"mean_coherence phil={m_phil:.3f} annotated={m_ann:.3f}")

    # sort plot by phylogeny order in SPECIES_META, keeping only discovered species
    order = [a for a in SPECIES_META if a in results]
    order += [a for a in results if a not in order]  # any unexpected extras last

    plot(results, order, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "functional_coherence.png"))

    out = {
        "params": {"min_cluster": MIN_CLUSTER, "recipe_threshold": RECIPE_THRESHOLD,
                   "metric": "mean pairwise GO Jaccard per cluster (PHILHARMONIC nb/02)"},
        "per_species": {
            a: {
                "label": SPECIES_META.get(a, (a, ""))[0],
                "phylum": SPECIES_META.get(a, ("", ""))[1],
                "n_clusters_total": r["n_clusters_total"],
                "n_scored_phil": r["n_scored_phil"],
                "n_clustered_proteins": r["n_clustered_proteins"],
                "annotation_coverage": r["annotation_coverage"],
                "mean_coherence_phil": float(np.mean(r["scores_phil"])) if r["scores_phil"] else None,
                "mean_coherence_annotated": float(np.mean(r["scores_annotated"])) if r["scores_annotated"] else None,
                "median_coherence_phil": float(np.median(r["scores_phil"])) if r["scores_phil"] else None,
                "scores_phil": r["scores_phil"],
                "scores_annotated": r["scores_annotated"],
            }
            for a, r in results.items()
        },
    }
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "functional_coherence.json")
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  wrote {outpath}")


if __name__ == "__main__":
    main()
