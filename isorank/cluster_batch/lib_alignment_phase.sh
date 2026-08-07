#!/bin/bash
# Shared body of the all-pairs IsoRank alignment phase -- sourced by both
# 05_run_all_pairs_alignment.sbatch (standalone) and 00_run_all.sbatch
# (orchestrator).
#
# Assumes: cwd is isorank/, venv activated, JOBS set. Requires
# 04_prepare_networks_and_pairs.sh to have already produced every
# net/<short>.tsv and net/<a>-<b>.tsv.

run_alignment_phase() {
    mkdir -p cluster_batch/logs

    mapfile -t SHORTS < <(tail -n +2 species.txt | cut -f1)
    : > cluster_batch/align_pairs.txt
    for ((i = 0; i < ${#SHORTS[@]}; i++)); do
        for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
            echo "${SHORTS[i]} ${SHORTS[j]}" >> cluster_batch/align_pairs.txt
        done
    done
    echo "found $(wc -l < cluster_batch/align_pairs.txt) pairs"

    run_one_alignment() {
        a="$1"; b="$2"
        if [[ -s "results/${a}_${b}_alignment_scores.json" ]]; then
            echo "skip: $a-$b"
            return 0
        fi
        if python3 3_isorank.py "$a" "$b" > "cluster_batch/logs/${a}_${b}.log" 2>&1; then
            echo "done:  $a-$b"
        else
            echo "FAILED: $a-$b (see cluster_batch/logs/${a}_${b}.log)"
        fi
    }
    export -f run_one_alignment

    xargs -P "$JOBS" -n 2 bash -c 'run_one_alignment "$@"' _ < cluster_batch/align_pairs.txt
}
