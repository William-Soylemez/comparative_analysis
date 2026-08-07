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
# LCC-restriction + bitscore similarity (2_similarity.py) -> all-pairs IsoRank
# alignment (3_isorank.py) -> summary table. Same fan-out logic as the
# standalone 01-06 scripts (shared via lib_blast_phase.sh /
# lib_alignment_phase.sh), just run back-to-back in one allocation instead of
# separate submissions, so JOBS = full node core count is available to every
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
# Vista has no "blast" module (confirmed via `module spider blast` -> not
# found) -- using NCBI's static BLAST+ binaries instead, unpacked once under
# vista/ncbi-blast-current (see cluster_batch README/chat history for the
# download+symlink one-liner). Adjust this path if you installed elsewhere.
export PATH="/work/11301/wsoylemez/vista/ncbi-blast-current/bin:$PATH"
source /work/11301/wsoylemez/vista/venv/bin/activate
JOBS="${SLURM_CPUS_PER_TASK:-144}"

echo "############ [1/5] proteomes ############"
bash 01_download_proteomes.sh "$FASTA_BASE_DIR"

echo "############ [2/5] blast DBs ############"
bash 02_build_blastdbs.sh

echo "############ [3/5] all-pairs reciprocal BLAST ############"
source lib_blast_phase.sh
run_blast_phase

echo "############ [4/5] LCC restriction + bitscore similarity ############"
bash 04_prepare_networks_and_pairs.sh

echo "############ [5/5] all-pairs IsoRank alignment + summary ############"
cd ..   # isorank/ -- 3_isorank.py reads net/... relative to here
# Re-activate the isorank venv (netalign/numpy/pandas live here, not in the
# general vista/venv sourced at the top of this script). Without this, the
# alignment phase would keep running in the general venv and fail with
# ModuleNotFoundError: netalign.
source ../venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source cluster_batch/lib_alignment_phase.sh
run_alignment_phase

python3 cluster_batch/06_aggregate_summary.py

echo "all done."
