#!/bin/bash
#SBATCH --job-name=isorank_all_pairs
#SBATCH --output=isorank_all_pairs_%j.out
#SBATCH --error=isorank_all_pairs_%j.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# One-job orchestrator running the entire all-pairs IsoRank pipeline
# end to end: proteomes -> blastdbs -> all-pairs reciprocal BLAST ->
# LCC-restriction + bitscore similarity -> all-pairs IsoRank alignment ->
# summary table -> decay-curve plots. Same fan-out logic as the standalone
# 01-07 scripts (shared via lib_blast_phase.sh / lib_alignment_phase.sh),
# just run back-to-back in one allocation instead of six separate
# submissions, so JOBS = full node core count is available to every
# parallel phase in turn -- phases are sequential (BLAST must finish before
# the alignment inputs can be built), but each phase individually still
# gets to fan out across all 144 cores via xargs -P, so there's no loss of
# parallelism vs. running the scripts one at a time yourself.
#
# Every phase is resume-safe (skips work already done), so if this times
# out or a step fails partway, just fix the issue and resubmit -- it picks
# up where it left off rather than redoing everything.
#
# Wall time: dominated by the BLAST phase (Rhizophagus's 9 pairs are ~90%
# of total BLAST search-space across all 45 pairs); --time=20:00:00 gives
# headroom over the earlier 12h/4h split-job estimates, but this is a new,
# untested scale -- if it times out, just resubmit this same script.
#
# Usage:
#   sbatch 00_run_all.sbatch [FASTA_BASE_DIR]
#
# FASTA_BASE_DIR (optional): a directory containing
#   <acc>_results/<acc>_unfiltered.fasta for each of the 6 new species'
# accessions -- copied straight in as proteomes/<short>.faa (guaranteed to
# match the IDs already in each species' network file). Falls back to an
# NCBI `datasets` download per-species if a given accession's fasta isn't
# found there (needs `datasets` on PATH and outbound internet -- most TACC
# compute nodes don't have this, so if you're relying on the fallback,
# run 01_download_proteomes.sh once on the login node first instead).

set -euo pipefail
FASTA_BASE_DIR="${1:-}"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

source /work/11301/wsoylemez/vista/philharmonic/sbatch_jobs/common.sh
module load gcc cuda python3
module load blast   # adjust if Vista uses a different module name -- check `module spider blast`
source /work/11301/wsoylemez/vista/venv/bin/activate
JOBS="${SLURM_CPUS_PER_TASK:-144}"

echo "############ [1/6] proteomes ############"
bash 01_download_proteomes.sh "$FASTA_BASE_DIR"

echo "############ [2/6] blast DBs ############"
bash 02_build_blastdbs.sh

echo "############ [3/6] all-pairs reciprocal BLAST ############"
source lib_blast_phase.sh
run_blast_phase

echo "############ [4/6] LCC restriction + bitscore similarity ############"
bash 04_prepare_networks_and_pairs.sh

echo "############ [5/6] all-pairs IsoRank alignment ############"
cd ..   # isorank/ -- run_full_alignment.py reads net/... relative to here
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source cluster_batch/lib_alignment_phase.sh
run_alignment_phase

echo "############ [6/6] summary + plots ############"
python3 cluster_batch/06_aggregate_summary.py
bash cluster_batch/07_plot_all_pairs.sh

echo "all done."
