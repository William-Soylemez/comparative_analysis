#!/bin/bash
#SBATCH --job-name=isorank_blast
#SBATCH --output=isorank_blast_%j.out
#SBATCH --error=isorank_blast_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# All-pairs reciprocal BLASTP across every species in species.tsv (10
# species -> C(10,2) = 45 pairs -> 90 directional blastp runs), fanned out
# across this node's cores with xargs -P, one blastp process per direction.
# Resume-safe: any blast_out/<a>_vs_<b>.tsv that already exists and is
# non-empty is skipped, so this can just be resubmitted if it times out or
# a handful of jobs fail.
#
# blastp flags (-evalue 1e-5, -outfmt 6 qseqid/sseqid/bitscore/evalue/pident,
# no -max_target_seqs cap) are reconstructed to match the columns
# build_rblast.py / build_rblast_bitscore.py already expect and the
# multi-hit-per-query pattern seen in the 4 species done earlier (e.g.
# ncrassa->calbicans: 33640 hits over 9757 queries) -- if you know the exact
# flags used for the original 4 pairs, swap them in here for consistency.
#
# Rhizophagus (30,147 proteins) dominates total runtime -- its 9 pairs are
# roughly 90% of the total BLAST search-space across all 45 pairs, so expect
# this job's wall time to be dominated by whichever of those 9 finish last.
#
# Usage:
#   sbatch 03_run_all_pairs_blast.sbatch

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

source /work/11301/wsoylemez/vista/philharmonic/sbatch_jobs/common.sh
module load gcc cuda python3
module load blast   # adjust if Vista uses a different module name -- check `module spider blast`
source /work/11301/wsoylemez/vista/venv/bin/activate
JOBS="${SLURM_CPUS_PER_TASK:-144}"

mkdir -p ../blast_out logs

# Build the full pair list (45 unordered pairs, both directions = 90 jobs).
: > pair_jobs.txt
mapfile -t SHORTS < <(tail -n +2 species.tsv | cut -f1)
for ((i = 0; i < ${#SHORTS[@]}; i++)); do
    for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
        echo "${SHORTS[i]} ${SHORTS[j]}" >> pair_jobs.txt
        echo "${SHORTS[j]} ${SHORTS[i]}" >> pair_jobs.txt
    done
done
N_PAIRS=$(( ${#SHORTS[@]} * (${#SHORTS[@]} - 1) / 2 ))
echo "found ${#SHORTS[@]} species -> $N_PAIRS pairs -> $(wc -l < pair_jobs.txt) directional blast jobs"

run_one() {
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
export -f run_one

xargs -P "$JOBS" -n 2 bash -c 'run_one "$@"' _ < pair_jobs.txt
