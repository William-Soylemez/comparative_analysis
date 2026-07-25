#!/bin/bash
#SBATCH --job-name=isorank_align
#SBATCH --output=isorank_align_%j.out
#SBATCH --error=isorank_align_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# Run IsoRank alignment (run_full_alignment.py) for all 45 species pairs,
# LCC-restricted networks + bitscore similarity (per the approved plan),
# fanned out across this node's cores with xargs -P. Each pair is an
# independent netalign compute_isorank() call.
#
# Resume-safe: skips any pair whose <a>_lcc_<b>_lcc_alignment_scores.json
# already exists, so this can just be resubmitted after a partial run.
#
# Memory: compute_isorank builds an n1xn1 and an n2xn2 transition matrix per
# pair. Actual PPI-network node counts (not proteome sizes) top out at
# rhizophagus's 18,537 -- its 9 pairs are the memory-heaviest (~2.75GB just
# for its own P1 matrix, ~7-8GB peak per pair including everything else);
# the other 36 pairs are under ~1.5GB each. Worst case, all 45 running at
# once is roughly ~100-120GB total, which should comfortably fit under
# --mem=0 on a full node -- but this is untested at this scale, so keep an
# eye on the first run.
#
# Usage:
#   sbatch 05_run_all_pairs_alignment.sbatch
# Requires 04_prepare_networks_and_pairs.sh to have already produced every
# net/<short>_lcc.tsv and net/<a>_lcc-<b>_lcc.tsv.

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR/.."   # isorank/ -- run_full_alignment.py reads net/... relative to here

source /work/11301/wsoylemez/vista/philharmonic/sbatch_jobs/common.sh
module load gcc cuda python3
source venv/bin/activate
# Force numpy/BLAS to single-threaded so xargs -P can use every core itself
# instead of each of the 45 alignment processes fighting over all of them.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-144}"

mkdir -p cluster_batch/logs

mapfile -t SHORTS < <(tail -n +2 cluster_batch/species.tsv | cut -f1)
: > cluster_batch/align_pairs.txt
for ((i = 0; i < ${#SHORTS[@]}; i++)); do
    for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
        echo "${SHORTS[i]} ${SHORTS[j]}" >> cluster_batch/align_pairs.txt
    done
done
echo "found $(wc -l < cluster_batch/align_pairs.txt) pairs"

run_one() {
    a="$1"; b="$2"
    prefix="${a}_lcc_${b}_lcc"
    if [[ -s "${prefix}_alignment_scores.json" ]]; then
        echo "skip: $a-$b"
        return 0
    fi
    if python3 run_full_alignment.py "${a}_lcc" "${b}_lcc" \
            > "cluster_batch/logs/${a}_${b}.log" 2>&1; then
        echo "done:  $a-$b"
    else
        echo "FAILED: $a-$b (see cluster_batch/logs/${a}_${b}.log)"
    fi
}
export -f run_one

xargs -P "$JOBS" -n 2 bash -c 'run_one "$@"' _ < cluster_batch/align_pairs.txt

python3 cluster_batch/06_aggregate_summary.py
