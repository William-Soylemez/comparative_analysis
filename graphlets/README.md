# graphlets

Count the six connected **4-node graphlets** — P4 (path), claw (3-star), C4
(cycle), paw (triangle+tail), diamond (K4−e), K4 (clique) — in each species' PPI
network, so their structural *composition* can be compared across the 10 fungal
species. The headline output is a two-panel figure: absolute counts and each
species' graphlet mix.

Counting uses **ORCA** (fast combinatorial orbit counting), which gives the same
exact counts as walking every 4-subgraph, without actually enumerating them.

## Pipeline

Run from this directory with the project venv (`../venv/bin/python`); species are
read from `../input/`, outputs go to `results/`.

1. **Build the ORCA binary once:** `cd orca_build && make` compiles `orca` from
   `orca.cpp` (vendored from the ORCA repo — see Attribution below). The binary
   is platform-specific and gitignored, so rebuild it after checkout.
2. **`orca_count4.py --prep`** reads each `../input/<acc>/<acc>_network.positive.tsv`,
   remaps proteins to `0..n-1`, dedupes edges, and writes an ORCA input file per
   species into `orca_build/work/` (plus a manifest).
3. **Run ORCA on each input** to produce per-node orbit counts:
   `for f in orca_build/work/*.in; do orca_build/orca node 4 "$f" "${f%.in}.orbits"; done`
4. **`orca_count4.py --parse`** sums the per-node orbits into global graphlet
   counts and writes `results/graphlet4_counts.json`.
5. **`plot_graphlet4.py`** draws the profile figure for three representative
   species (sparse → dense, different phyla).

## Inputs

- `../input/<acc>/<acc>_network.positive.tsv` — PPI edge list (protein-name pairs;
  weights ignored — graphlets are purely topological).

## Outputs (in `results/`)

- `graphlet4_counts.json` — per-species node/edge counts and the six 4-graphlet
  totals (plus P3 path and triangle for reference).
- `graphlet4_profiles.png` — absolute counts (log) + composition (% of each
  species' 4-graphlets) for the three exemplars.

## Scripts & files

- `orca_count4.py` — the `--prep` / `--parse` driver (steps 2 and 4).
- `plot_graphlet4.py` — the profile figure; its 3 exemplar species are hardcoded.
- `orca_build/` — `orca.cpp` (source) + `Makefile`; the compiled `orca` and the
  `work/` intermediates are gitignored.

## Design decisions

- **ORCA, not exhaustive enumeration.** Orbit counting is combinatorial (orders
  of magnitude faster) yet gives exact counts — enumerating every induced
  4-subgraph is infeasible on the denser species (millions of edges).
- **Orbit → graphlet by node-multiplicity.** ORCA reports per-node *orbit*
  counts; each graphlet total is `sum(orbit) / (nodes of that orbit per
  graphlet)` — e.g. C4 = `sum(orbit8)/4`, claw = `sum(orbit7)/1`.
- **Unweighted / topological.** D-SCRIPT edge weights are dropped; graphlets are
  about wiring shape, not edge strength.
- **Three exemplars for the figure**, not all 10 — a readable sparse→mid→dense
  comparison across three phyla, rather than a crowded 10-series plot.

## Attribution

`orca.cpp` is vendored from ORCA (ORbit Counting Algorithm). To rebuild from
upstream, clone the repo and use its `orca.cpp`:

- Repo: https://github.com/thocevar/orca
- Paper: Hočevar T. & Demšar J., *A combinatorial approach to graphlet counting*,
  Bioinformatics 30(4):559–565 (2014). doi:10.1093/bioinformatics/btu245

