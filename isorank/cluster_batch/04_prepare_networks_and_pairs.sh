#!/bin/bash
# Two quick prep steps once BLAST is done, before the alignment batch:
#   1. Extract net/<short>.tsv (2-column edge list) from each species'
#      <acc>_network.positive.tsv (which has a 3rd D-SCRIPT-score column
#      that compute_adjacency() can't handle -- see run_full_alignment.py).
#   2. For every one of the 45 pairs: restrict both species' networks to
#      their largest connected component (make_lcc_networks.py, writing
#      net/<short>_lcc.tsv) and build the bitscore-based similarity file
#      (build_rblast_bitscore.py, writing net/<a>_lcc-<b>_lcc.tsv) -- so
#      run_full_alignment.py can be pointed at "<short>_lcc" for every
#      species with no further changes, exactly as done for the first 3
#      pairs (LCC-restricted, bitscore-scored, per the approved plan).
#
# Cheap (seconds), no sbatch needed -- run directly on the login node or as
# a quick srun/interactive step.
#
# Usage: bash 04_prepare_networks_and_pairs.sh
# Expects species.tsv, ../../<species_dir>/<acc>_network.positive.tsv,
#         ../blast_out/<a>_vs_<b>.tsv (from 03_run_all_pairs_blast.sbatch)
# Writes ../net/<short>.tsv, ../net/<short>_lcc.tsv, ../net/<a>_lcc-<b>_lcc.tsv

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p ../net
BASE_DIR=../..  # comparative_analysis/, containing <species_dir>/<acc>_network.positive.tsv

echo "== extracting 2-column edge lists =="
tail -n +2 species.tsv | while IFS=$'\t' read -r short acc organism species_dir; do
    out="../net/${short}.tsv"
    if [[ -s "$out" ]]; then
        echo "skip: $short -- $out already exists"
        continue
    fi
    src="$BASE_DIR/$species_dir/${acc}_network.positive.tsv"
    if [[ ! -s "$src" ]]; then
        echo "  ! missing $src, skipping $short" >&2
        continue
    fi
    cut -f1,2 "$src" > "$out"
    echo "wrote $out ($(wc -l < "$out") edges)"
done

echo "== restricting every species to its largest connected component =="
mapfile -t SHORTS < <(tail -n +2 species.tsv | cut -f1)
( cd .. && source venv/bin/activate && python3 make_lcc_networks.py "${SHORTS[@]}" )

echo "== building bitscore-based similarity files for all 45 pairs =="
( cd .. && source venv/bin/activate
  for ((i = 0; i < ${#SHORTS[@]}; i++)); do
      for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
          a="${SHORTS[i]}"; b="${SHORTS[j]}"
          out="net/${a}_lcc-${b}_lcc.tsv"
          if [[ -s "$out" ]]; then
              echo "skip: $a-$b -- $out already exists"
              continue
          fi
          python3 build_rblast_bitscore.py "$a" "$b" "${a}_lcc" "${b}_lcc"
      done
  done
)
