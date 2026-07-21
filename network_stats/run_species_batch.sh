#!/bin/bash
#SBATCH --job-name=network_stats
#SBATCH --output=network_stats_%j.out
#SBATCH --error=network_stats_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# Compute global network stats for every species under BASE_DIR, one
# single-threaded compute_network_stats.py process per species, fanned out
# across this node's cores with xargs -P. Each species job is independent
# (unlike the collaborator's post_cluster_par.sh, there's no per-species
# internal parallelism to hand off to snakemake here -- this is just
# embarrassingly parallel across species).
#
# Usage:
#   sbatch run_species_batch.sh <species_base_dir> [out_dir] [accessions_file]
#
# With no accessions_file: every immediate subdirectory of <species_base_dir>
# that has a *_network.positive.tsv file in it is treated as a species dir.
#
# With accessions_file (one accession per line, blank lines and #comments
# ignored): only <species_base_dir>/<acc>_results is used for each accession
# -- matching the philharmonic pipeline's <acc>_results naming convention.

set -euo pipefail

BASE_DIR="${1:?Usage: sbatch run_species_batch.sh <species_base_dir> [out_dir] [accessions_file]}"
OUT_DIR="${2:-$BASE_DIR/network_stats_out}"
ACCESSIONS_FILE="${3:-}"
# sbatch copies this script into /var/spool/slurmd/job<id>/ and runs it from
# there, so BASH_SOURCE's own dirname is useless under Slurm -- use
# SLURM_SUBMIT_DIR (cwd at submission time) instead, which is where
# compute_network_stats.py / aggregate_stats.py actually live.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
JOBS="${SLURM_CPUS_PER_TASK:-144}"

mkdir -p "$OUT_DIR/logs"

if [[ -n "$ACCESSIONS_FILE" ]]; then
    # One species dir per accession: <base_dir>/<acc>_results, skipping any
    # that don't exist or lack a network file (with a warning).
    : > "$OUT_DIR/species_list.txt"
    while read -r acc; do
        [[ -z "$acc" || "$acc" == \#* ]] && continue
        d="$BASE_DIR/${acc}_results"
        if compgen -G "$d/*_network.positive.tsv" > /dev/null; then
            echo "$d" >> "$OUT_DIR/species_list.txt"
        else
            echo "  ! skipping $acc: no *_network.positive.tsv in $d" >&2
        fi
    done < "$ACCESSIONS_FILE"
else
    # One species dir per line: any immediate subdirectory of BASE_DIR that
    # has a *_network.positive.tsv file in it.
    find "$BASE_DIR" -maxdepth 1 -mindepth 1 -type d \
        | while read -r d; do
            compgen -G "$d/*_network.positive.tsv" > /dev/null && echo "$d"
          done \
        > "$OUT_DIR/species_list.txt"
fi

N_SPECIES=$(wc -l < "$OUT_DIR/species_list.txt")
echo "found $N_SPECIES species dirs under $BASE_DIR, running with -P $JOBS"

run_one() {
    dir="$1"
    acc="$(basename "$dir")"
    if python3 "$SCRIPT_DIR/compute_network_stats.py" "$dir" --quiet \
            --out "$OUT_DIR/${acc}_stats.json" \
            > "$OUT_DIR/logs/${acc}.log" 2>&1; then
        echo "done:  $acc"
    else
        echo "FAILED: $acc (see $OUT_DIR/logs/${acc}.log)"
    fi
}
export -f run_one
export SCRIPT_DIR OUT_DIR

xargs -I{} -P "$JOBS" bash -c 'run_one "$@"' _ {} < "$OUT_DIR/species_list.txt"

python3 "$SCRIPT_DIR/aggregate_stats.py" "$OUT_DIR"
