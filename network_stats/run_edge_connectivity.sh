#!/bin/bash
#SBATCH --job-name=cluster_edge_conn
#SBATCH --output=cluster_edge_conn_%j.out
#SBATCH --error=cluster_edge_conn_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# Compute cluster_graph.avg_edge_connectivity_largest_cc for every species
# under BASE_DIR that doesn't already have it (a brand-new field -- missing
# for all 79 current species; networkx has no built-in
# average_edge_connectivity, so this uses the same parallel machinery as
# avg_node_connectivity, swapping in local_edge_connectivity +
# build_auxiliary_edge_connectivity -- see compute_cluster_graph_stats.py).
# Same asymptotic cost as node connectivity (~similar CPU-hours), though
# edge connectivity's max-flow skips node connectivity's node-splitting
# trick so it may come in cheaper in practice. --skip avg_node_connectivity
# keeps this job from also (re)computing node connectivity -- that's
# run_node_connectivity.sh's job, so the two can be scheduled/rerun
# completely independently of each other.
#
# Unlike run_species_batch.sh, this deliberately does NOT fan out with
# `xargs -P` (one species per core) -- that model is right for uniform,
# cheap, per-species stats, but wrong here: the work per species is wildly
# uneven (minutes to potentially 90+ hours), so fanning across species would
# leave most cores idle while the single biggest species finishes alone on
# one core. Instead this gives EVERY core to ONE species at a time (--procs
# = the whole node), processing the most expensive species first.
#
# Safe to run independently of run_species_batch.sh / run_node_connectivity.sh,
# in any order, any number of times (e.g. after adding new species): all
# three write into the same --out files, each skips whatever's already
# there, and compute_network_stats.py re-reads + merges right before writing
# so a concurrently-running job computing a different field can't have its
# result clobbered (see merge_results in compute_network_stats.py). Safe to
# requeue after a timeout too -- already-finished species are an instant
# no-op and only the remaining queue gets worked through.
#
# Usage:
#   sbatch run_edge_connectivity.sh <species_base_dir> <out_dir> [accessions_file]
#
# <out_dir> should be the SAME directory used with run_species_batch.sh (it
# will be created if this is the very first job run against a fresh
# species_base_dir).

set -euo pipefail

BASE_DIR="${1:?Usage: sbatch run_edge_connectivity.sh <species_base_dir> <out_dir> [accessions_file]}"
OUT_DIR="${2:?Usage: sbatch run_edge_connectivity.sh <species_base_dir> <out_dir> [accessions_file]}"
ACCESSIONS_FILE="${3:-}"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source /work/11301/wsoylemez/vista/philharmonic/sbatch_jobs/common.sh
module load gcc cuda python3
source /work/11301/wsoylemez/vista/venv/bin/activate
JOBS="${SLURM_CPUS_PER_TASK:-144}"

mkdir -p "$OUT_DIR/logs"

# Print exactly what was resolved before doing anything else -- a bad path,
# an empty/missing accessions file, or a directory-naming mismatch should be
# obvious from this instead of inferred later from "found 0 species dirs".
echo "resolved cwd: $(pwd)"
echo "BASE_DIR: $BASE_DIR"
echo "OUT_DIR: $OUT_DIR"
if [[ -n "$ACCESSIONS_FILE" ]]; then
    if [[ -f "$ACCESSIONS_FILE" ]]; then
        echo "ACCESSIONS_FILE: $ACCESSIONS_FILE ($(wc -l < "$ACCESSIONS_FILE") lines)"
    else
        echo "ACCESSIONS_FILE: $ACCESSIONS_FILE -- DOES NOT EXIST at this path"
    fi
else
    echo "ACCESSIONS_FILE: <none, auto-discovering all species under BASE_DIR>"
fi

if [[ -n "$ACCESSIONS_FILE" ]]; then
    : > "$OUT_DIR/edge_conn_species_list.txt"
    while read -r acc; do
        acc="${acc%$'\r'}"  # strip a trailing \r in case the file has Windows line endings
        [[ -z "$acc" || "$acc" == \#* ]] && continue
        found=""
        for d in "$BASE_DIR/${acc}_results" "$BASE_DIR/${acc}"; do
            if compgen -G "$d/*_network.positive.tsv" > /dev/null; then
                echo "$d" >> "$OUT_DIR/edge_conn_species_list.txt"
                found=1
                break
            fi
        done
        [[ -z "$found" ]] && echo "  ! skipping $acc: no *_network.positive.tsv in $BASE_DIR/${acc}_results or $BASE_DIR/${acc}" >&2
    done < "$ACCESSIONS_FILE"
else
    find "$BASE_DIR" -maxdepth 1 -mindepth 1 -type d \
        | while read -r d; do
            compgen -G "$d/*_network.positive.tsv" > /dev/null && echo "$d"
          done \
        > "$OUT_DIR/edge_conn_species_list.txt"
fi

# Sort species dirs by descending recorded pair count (from a prior run's
# meta -- either connectivity job's meta works, since total_pairs is the
# same C(largest_cc_size, 2) regardless of which metric it's attached to) so
# the most expensive species run first -- if the job times out, only cheap
# species are left for a quick follow-up run. Species never touched before
# (no recorded meta at all) sort last, since their cost is unknown.
python3 - "$OUT_DIR" < "$OUT_DIR/edge_conn_species_list.txt" > "$OUT_DIR/edge_conn_species_list.sorted.txt" <<'PY'
import json, os, sys

out_dir = sys.argv[1]
rows = []
for line in sys.stdin:
    d = line.strip()
    if not d:
        continue
    acc = os.path.basename(d)
    stats_path = os.path.join(out_dir, f"{acc}_stats.json")
    pairs = 0
    try:
        with open(stats_path) as f:
            cg = json.load(f).get("cluster_graph") or {}
        meta = cg.get("avg_edge_connectivity_meta") or cg.get("avg_node_connectivity_meta") or {}
        pairs = meta.get("total_pairs", 0) or 0
    except (OSError, json.JSONDecodeError):
        pass
    rows.append((pairs, d))

rows.sort(key=lambda r: r[0], reverse=True)
for _, d in rows:
    print(d)
PY
mv "$OUT_DIR/edge_conn_species_list.sorted.txt" "$OUT_DIR/edge_conn_species_list.txt"

N_SPECIES=$(wc -l < "$OUT_DIR/edge_conn_species_list.txt")
echo "found $N_SPECIES species dirs under $BASE_DIR, running sequentially with --procs $JOBS each"

t_start=$(date +%s)
i=0
while read -r dir; do
    i=$((i + 1))
    acc="$(basename "$dir")"
    out="$OUT_DIR/${acc}_stats.json"
    echo "[$i/$N_SPECIES] $acc ..."
    if python3 "$SCRIPT_DIR/compute_network_stats.py" "$dir" --quiet \
            --skip avg_node_connectivity \
            --out "$out" --procs "$JOBS" \
            > "$OUT_DIR/logs/${acc}_edge_conn.log" 2>&1; then
        echo "  done:  $acc"
    else
        echo "  FAILED: $acc (see $OUT_DIR/logs/${acc}_edge_conn.log)"
    fi
    echo "  elapsed so far: $(( $(date +%s) - t_start ))s"
done < "$OUT_DIR/edge_conn_species_list.txt"

python3 "$SCRIPT_DIR/aggregate_stats.py" "$OUT_DIR"
