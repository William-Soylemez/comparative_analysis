# isorank

Pairwise **network alignment** between two species' PPI networks with IsoRank:
BLAST the proteomes against each other, turn the reciprocal hits into a
sequence-similarity signal (normalized bitscore), and let IsoRank blend that
with network topology to match proteins across species. Every matched pair keeps
its IsoRank score, so a confidence cutoff can be chosen downstream.

## Pipeline

Three scripts, run **from this directory** for a pair of species (short names,
e.g. `scer calbicans`) with the project venv (`../venv/bin/python`):

1. **`1_blast.sh <a> <b>`** — build a protein BLAST DB per species (from
   `proteomes/<short>.faa`) and run reciprocal `blastp` both directions →
   `blast_out/<a>_vs_<b>.tsv`, `<b>_vs_<a>.tsv`.
2. **`2_similarity.py <a> <b>`** — restrict each network to its largest connected
   component → `net/<short>.tsv`; merge both BLAST directions (max bitscore per
   protein pair), min-max normalize to [0, 1] → `net/<a>-<b>.tsv`.
3. **`3_isorank.py <a> <b>`** — run IsoRank (via the netalign library) on the two
   networks + the similarity signal, greedily match **100% of the smaller
   network's nodes**, and write every pair with its score →
   `results/<a>_<b>_alignment_scored.tsv` + `_scores.json`.

`species.txt` maps each short name to its NCBI accession, so the scripts can find
that species' proteome and network.

## Inputs

- `proteomes/<short>.faa` — protein FASTA (the exact set PHILHARMONIC used, so
  IDs match the network). Acquired separately (see `cluster_batch/01`).
- `../input/<acc>/<acc>_network.positive.tsv` — the species' PPI network.

## Outputs (in `results/`)

- `<a>_<b>_alignment_scored.tsv` — one row per node of the smaller network:
  `a-protein  b-protein  score`.
- `<a>_<b>_alignment_scores.json` — summary (pair count, coverage, score
  distribution, IsoRank params).

## Design decisions

- **Normalized bitscore as the similarity metric**, not `-log10(evalue)`. Both
  were tried; bitscore is the committed signal (max over both BLAST directions,
  then min-max normalized to [0, 1] within the pair).
- **Restrict to the largest connected component** before aligning. Most of these
  networks are already a single component; this just avoids IsoRank quirks on the
  occasional tiny disconnected island.
- **Match all of the smaller network, cut off later.** The greedy assignment
  covers 100% of the smaller species and keeps every score — a confidence
  threshold is a downstream choice, not baked into the alignment.
- **IsoRank via netalign's `compute_*` functions directly**, not the `netalign
  isorank` CLI, because the CLI returns only index pairs — we need the per-match
  score. Blend factor `ALPHA = 0.7` (topology-weighted), `NITER = 20`.

## Running all pairs on the cluster

`cluster_batch/` fans the same three scripts across all C(10,2) = 45 species
pairs on SLURM, resume-safe (each phase skips work already done):

- `01_download_proteomes.sh` → proteomes; `02_build_blastdbs.sh` + reciprocal
  BLAST (`03_*.sbatch` / `lib_blast_phase.sh`).
- `04_prepare_networks_and_pairs.sh` runs `2_similarity.py` for every pair;
  `05_*.sbatch` / `lib_alignment_phase.sh` runs `3_isorank.py` for every pair.
- `06_aggregate_summary.py` collects the per-pair `_scores.json` into
  `results/all_pairs_summary.{tsv,json}`.
- `00_run_all.sbatch` runs the whole thing end-to-end in one allocation.

It reads the same `species.txt` and produces the same per-pair results as running
the three scripts by hand — just across all 45 pairs at once. The `.sbatch` files
carry cluster-specific paths (TACC/Vista) you'll adjust for your environment.
