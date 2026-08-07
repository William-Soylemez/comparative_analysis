#!/bin/bash
#SBATCH --job-name=isorank_blast
#SBATCH --output=isorank_blast_%j.out
#SBATCH --error=isorank_blast_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=144
#SBATCH --mem=0

# All-pairs reciprocal BLASTP across every species in species.txt (10
# species -> C(10,2) = 45 pairs -> 90 directional blastp runs), fanned out
# across this node's cores with xargs -P. Fan-out logic lives in
# lib_blast_phase.sh (shared with 00_run_all.sbatch). Resume-safe: any
# blast_out/<a>_vs_<b>.tsv that already exists and is non-empty is skipped.
#
# blastp flags (-evalue 1e-5, -outfmt 6 qseqid/sseqid/bitscore/evalue/pident,
# no -max_target_seqs cap) match the columns 2_similarity.py expects and the
# multi-hit-per-query pattern seen in the 4 species done earlier (e.g.
# ncrassa->calbicans: 33640 hits over 9757 queries). Same flags as the
# single-pair 1_blast.sh.
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
# Vista has no "blast" module (confirmed via `module spider blast` -> not
# found) -- using NCBI's static BLAST+ binaries instead, unpacked once under
# vista/ncbi-blast-current (see cluster_batch README/chat history for the
# download+symlink one-liner). Adjust this path if you installed elsewhere.
export PATH="/work/11301/wsoylemez/vista/ncbi-blast-current/bin:$PATH"
source /work/11301/wsoylemez/vista/venv/bin/activate
JOBS="${SLURM_CPUS_PER_TASK:-144}"

source lib_blast_phase.sh
run_blast_phase
