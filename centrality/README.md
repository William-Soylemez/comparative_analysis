# centrality

Are high-**betweenness** ("connector") proteins — those on many shortest paths —
enriched for particular functions, consistently across the 10 fungal species?
The headline output is a GO-category × species heatmap of enrichment among
connector proteins.

## Pipeline

Run all scripts **from this directory** with the project venv
(`../venv/bin/python`). They discover species under `../input/`, read the
ontology from `../go_data/`, and write everything to `results/`.

1. **`run_all_betweenness.py`** computes betweenness for every node in every
   species and overlays the per-species "decay curves" (betweenness vs. top
   percentile). Eyeball the overlay to confirm the elbow is consistent — this is
   what justifies a single **top-5%** connector cutoff.
   (`betweenness_threshold_plot.py` is the single-species version + shared library.)
2. **`fisher_enrichment.py`** takes the top 5% of nodes as set *H*, and for each
   GO term runs a one-sided Fisher exact test for over-representation in *H* vs
   the rest of the annotated network, BH-FDR corrected. Add `--slim` (GO-slim
   roll-up), `--propagate` (ancestors), or `--centrality degree` as variants.
3. **`plot_enrichment_heatmap.py <tag>`** draws the heatmap for one variant,
   showing fold enrichment only where significant (`q ≤ 0.05`).

```bash
../venv/bin/python run_all_betweenness.py
../venv/bin/python fisher_enrichment.py --slim
../venv/bin/python plot_enrichment_heatmap.py slim
```

## Inputs

- `../input/<acc>/<acc>_network.positive.tsv` — PPI edge list; weights are
  D-SCRIPT similarity scores.
- `../input/<acc>/<acc>_GO_map.csv` — per-protein GO annotations.
- `../go_data/{go-basic,goslim_generic}.obo` + `../goslim_util.py` — ontology and
  slim roll-up.

## Outputs (in `results/`, by `<tag>` = `direct` | `slim` | `propagated`, or `degree_*`)

- `<acc>_betweenness_decay.{png,json}`, `all_species_betweenness_decay.png`,
  `all_species_betweenness_summary.json` — threshold analysis.
- `<acc>_enrichment_<tag>.tsv`, `enrichment_summary_<tag>.json` — per-term stats.
- `enrichment_<tag>_heatmap.png` — the heatmap.
- `<acc>_node_{betweenness,degree}.json` — score **caches**; deletable, regenerated on next run.

## Design decisions

- **Unweighted betweenness.** D-SCRIPT weights are *similarities*, but graph libs
  treat weights as *distances* — using them raw would invert the geometry. Scores
  are size-normalized so cutoffs compare across species.
- **One top-5% cutoff for all species** (justified by the decay-curve overlay),
  not a per-species elbow.
- ***H* is defined over all nodes, then intersected with the annotated set** — a
  purely topological "connector" definition; the annotated shortfall is reported,
  not hidden.
- **GO-slim is the lens for cross-species claims.** Direct/propagated terms carry
  a real confound (metazoan terms leaking into fungal annotations, uneven
  annotation depth). Trust `enrichment_slim_heatmap.png` for the headline; treat
  `direct` as a drill-down. `--centrality degree` is a control for whether signal
  is specific to brokerage vs. plain hubness.
