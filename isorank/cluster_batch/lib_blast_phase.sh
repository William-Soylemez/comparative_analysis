#!/bin/bash
# Shared body of the all-pairs reciprocal BLASTP phase -- sourced by both
# 03_run_all_pairs_blast.sbatch (standalone) and 00_run_all.sbatch
# (orchestrator) so the fan-out logic only lives in one place.
#
# Assumes: cwd is cluster_batch/, venv activated, blastp on PATH, JOBS set.

run_blast_phase() {
    mkdir -p ../blast_out logs

    : > pair_jobs.txt
    mapfile -t SHORTS < <(tail -n +2 species.txt | cut -f1)
    for ((i = 0; i < ${#SHORTS[@]}; i++)); do
        for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
            echo "${SHORTS[i]} ${SHORTS[j]}" >> pair_jobs.txt
            echo "${SHORTS[j]} ${SHORTS[i]}" >> pair_jobs.txt
        done
    done
    N_PAIRS=$(( ${#SHORTS[@]} * (${#SHORTS[@]} - 1) / 2 ))
    echo "found ${#SHORTS[@]} species -> $N_PAIRS pairs -> $(wc -l < pair_jobs.txt) directional blast jobs"

    run_one_blast() {
        q="$1"; s="$2"
        out="../blast_out/${q}_vs_${s}.tsv"
        if [[ -s "$out" ]]; then
            echo "skip: ${q}_vs_${s}"
            return 0
        fi
        if blastp -query "../proteomes/${q}.faa" -db "../blastdb/${s}" \
                -evalue 1e-5 -outfmt "6 qseqid sseqid bitscore evalue pident" \
                -out "$out.tmp" > "logs/${q}_vs_${s}.log" 2>&1; then
            mv "$out.tmp" "$out"
            echo "done:  ${q}_vs_${s} ($(wc -l < "$out") hits)"
        else
            echo "FAILED: ${q}_vs_${s} (see logs/${q}_vs_${s}.log)"
            rm -f "$out.tmp"
        fi
    }
    export -f run_one_blast

    xargs -P "$JOBS" -n 2 bash -c 'run_one_blast "$@"' _ < pair_jobs.txt
}
