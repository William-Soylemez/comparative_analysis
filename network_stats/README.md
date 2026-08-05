# network_stats

Global topology statistics for each species' PHILHARMONIC PPI network, computed
on **two graphs**: the protein graph (the network itself) and a derived
"cluster-of-clusters" graph. Per-species results are aggregated into one JSON
and rendered as an interactive comparison page.

## Pipeline

The compute scripts take a species dir as an argument and are designed to run on
a SLURM cluster over many species (the `run_*.sh` wrappers), so `BASE_DIR` /
`OUT_DIR` are parameters, not hardcoded. Locally, species live under `../input/`
and the viz build uses the project venv (`../venv/bin/python`).

1. **`compute_network_stats.py <species_dir> --out <acc>_stats.json`** — protein-graph
   stats, plus a nested `cluster_graph` block from `compute_cluster_graph_stats.py`.
   Field-level resume: rerunning only computes keys missing from an existing
   `--out`, so adding a new stat later backfills just that stat.
2. **`run_species_batch.sh`** — fans the per-species compute across cores, one
   core each (`xargs -P`). Every stat is cheap and single-threaded.
3. **`aggregate_stats.py <out_dir>`** — combines all `<acc>_stats.json` into
   `all_species_data.json` (atomic write; skips its own output file).
4. **`build_viz.py`** — embeds `all_species_data.json` into
   `species_viz.template.html` → `results/species_viz.html` (see below).

## Inputs

- `../input/<acc>/<acc>_network.positive.tsv` — PPI edge list.
- `../input/<acc>/<acc>_clusters.json` — PHILHARMONIC clusters (`members` lists).

## Outputs

- `<out_dir>/<acc>_stats.json` — per species; `all_species_data.json` — aggregated.
- `results/species_viz.html` — the interactive page (embedded data + live-fetch fallback).

## The stats

**Protein graph:** `node_count`, `edge_count`, `density`,
`num_connected_components`, `nodes_outside_largest_cc`, `largest_cc_size`,
`median_degree`, `num_philharmonic_clusters`, `global_clustering_coefficient`,
`num_bridges`. (Diameter and modularity are computed **only for the cluster
graph** — see below.)

**Cluster graph** (nodes = clusters; edge when ≥ `MIN_CROSSING_EDGES` = 50
protein–protein edges cross between two clusters' members; built from
`clusters.json`, *not* the PHILHARMONIC-provided `cluster_graph.tsv`, which used
pre-ReCIPE membership and didn't match): the same size/clustering/bridge stats,
plus `node_connectivity` / `edge_connectivity` (largest CC), `diameter`
(**exact** here — the graph is small), and `modularity` (the cluster graph's
**own** greedy communities).

## The viz

`species_viz.template.html` is the source (data replaced by an
`__EMBEDDED_DATA__` marker); `build_viz.py` splices `all_species_data.json` in.
The built page works two ways: double-clicked (`file://`) it renders the
embedded snapshot; served (`python3 -m http.server`) it also fetches
`all_species_data.json` from the same folder and prefers that live copy. Rebuild
after new stats with `../venv/bin/python build_viz.py`.

## Todo: average connectivity

Average connectivity was quite slow for the longest species and needs more aggressive
intra-species parallelization. I took a stab at this but it wasn't really working, so
a bit of tweaking will be needed.