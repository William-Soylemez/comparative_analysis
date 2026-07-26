"""
For GO:0016126 (ergosterol biosynthetic process), count per fungal species:
  - how many proteins are annotated with the term (from *_GO_map.csv)
  - how many distinct PHILHARMONIC clusters those proteins fall in (from *_clusters.json 'members')

Run from anywhere; paths are relative to the comparative_analysis root.
"""
import csv
import json
import sys
from pathlib import Path

GO_TERM = "GO:0016126"
ROOT = Path(__file__).resolve().parent.parent

# GCF_000146465.1 (E. intestinalis) excluded at user's request.
SPECIES = {
    "GCF_000146045.2": "Saccharomyces cerevisiae",
    "GCF_000149245.1": "Cryptococcus neoformans",
    "GCF_000203795.2": "Batrachochytrium dendrobatidis",
    "GCF_026210795.1": "Rhizophagus irregularis",
    "GCF_025024165.1": "Kickxella alabastrina",
    "GCA_014872705.1": "Agaricus bisporus",
    "GCF_000182965.3": "Candida albicans",
    "GCF_000182925.2": "Neurospora crassa",
    "GCA_025594325.1": "Blastocladiella emersonii",
    "GCF_025094135.1": "Mucor mucedo",
}


def annotated_proteins(acc):
    path = ROOT / acc / f"{acc}_GO_map.csv"
    hits = set()
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            go_list = row["GO_list"].split(";") if row["GO_list"] else []
            if GO_TERM in go_list:
                hits.add(row["prot_id"])
    return hits


def clusters_path(acc):
    final = ROOT / acc / f"{acc}_clusters.json"
    if final.exists():
        return final, False
    pre = ROOT / acc / f"{acc}_clusters.pre.json"
    if pre.exists():
        return pre, True
    raise FileNotFoundError(f"no clusters.json or clusters.pre.json for {acc}")


def clusters_containing(acc, proteins):
    path, is_pre = clusters_path(acc)
    clusters = json.load(open(path))
    hit_clusters = {}
    for cid, c in clusters.items():
        members_hit = proteins & set(c["members"])
        if members_hit:
            hit_clusters[cid] = sorted(members_hit)
    return hit_clusters, is_pre


def main():
    results = []
    for acc, organism in SPECIES.items():
        proteins = annotated_proteins(acc)
        if proteins:
            hit_clusters, is_pre = clusters_containing(acc, proteins)
        else:
            hit_clusters, is_pre = {}, False
        results.append(
            {
                "accession": acc,
                "organism": organism,
                "n_annotated_proteins": len(proteins),
                "n_clusters": len(hit_clusters),
                "clusters_pre_llm_naming": is_pre,
                "proteins": sorted(proteins),
                "cluster_membership": hit_clusters,
            }
        )

    out_path = ROOT / "ergosterol_biosynthesis" / "ergosterol_biosynthesis_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"{'Organism':<32}{'Accession':<18}{'GO-annotated proteins':<24}{'Distinct clusters':<20}{'Source'}")
    for r in results:
        src = "clusters.pre.json" if r["clusters_pre_llm_naming"] else "clusters.json"
        print(f"{r['organism']:<32}{r['accession']:<18}{r['n_annotated_proteins']:<24}{r['n_clusters']:<20}{src}")
    print(f"\nFull detail written to {out_path}")


if __name__ == "__main__":
    main()
