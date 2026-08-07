#!/bin/bash
# Build the IsoRank inputs (LCC-restricted networks + normalized-bitscore
# similarity) for all species pairs, by running 2_similarity.py once per pair.
# Cheap (seconds per pair) -- run after BLAST, before the alignment phase.
# Resume-safe: skips a pair whose net/<a>-<b>.tsv already exists.
#
# Usage: bash 04_prepare_networks_and_pairs.sh
# Expects ../species.txt, ../blast_out/<a>_vs_<b>.tsv (both directions, from
#         the BLAST phase), ../../input/<acc>/<acc>_network.positive.tsv
# Writes  ../net/<short>.tsv (LCC edge lists) and ../net/<a>-<b>.tsv (similarity)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."   # isorank/  (2_similarity.py + species.txt live here)
source ../venv/bin/activate

mapfile -t SHORTS < <(tail -n +2 species.txt | cut -f1)
for ((i = 0; i < ${#SHORTS[@]}; i++)); do
    for ((j = i + 1; j < ${#SHORTS[@]}; j++)); do
        a="${SHORTS[i]}"; b="${SHORTS[j]}"
        if [[ -s "net/${a}-${b}.tsv" ]]; then
            echo "skip: $a-$b -- net/${a}-${b}.tsv already exists"
            continue
        fi
        python3 2_similarity.py "$a" "$b"
    done
done
