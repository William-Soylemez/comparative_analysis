#!/bin/bash
# Generate a decay-curve plot (plot_score_histogram.py) for every one of
# the 45 species pairs aligned by 05_run_all_pairs_alignment.sbatch. Cheap
# (matplotlib, no numpy-heavy work) -- fine on the login node, or as the
# last step of 00_run_all.sbatch.
#
# Resume-safe in the trivial sense that it just regenerates every plot each
# time it's run (a few seconds each); doesn't skip existing ones since
# plots are cheap and you may want to rerun after a labeling/style tweak.
#
# Usage: bash 07_plot_all_pairs.sh
# Expects ../<a>_lcc_bs_<b>_lcc_bs_alignment_scored.tsv for each pair.
# Writes ../<a>_lcc_bs_<b>_lcc_bs_score_decay.png (+ .json) for each pair.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."   # isorank/
source ../venv/bin/activate

mapfile -t SHORTS < <(tail -n +2 cluster_batch/species.txt | cut -f1)
n_ok=0
n_missing=0
for ((i = 0; i < ${#SHORTS[@]}; i++)); do
    for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
        a="${SHORTS[i]}_lcc_bs"; b="${SHORTS[j]}_lcc_bs"
        if [[ ! -s "${a}_${b}_alignment_scored.tsv" ]]; then
            echo "  ! missing ${a}_${b}_alignment_scored.tsv, skipping" >&2
            n_missing=$((n_missing + 1))
            continue
        fi
        python3 plot_score_histogram.py "$a" "$b"
        n_ok=$((n_ok + 1))
    done
done
echo "plotted $n_ok pairs, $n_missing missing"
