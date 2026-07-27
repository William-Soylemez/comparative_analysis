#!/bin/bash
# Populate proteomes/<short>.faa for every species in species.tsv that
# doesn't already have one, preferring the exact protein set PHILHARMONIC
# itself used to build that species' network (so IDs are guaranteed to
# match net/<short>.tsv 1:1) over a fresh NCBI download.
#
# Usage: bash 01_download_proteomes.sh [FASTA_BASE_DIR]
#
# If FASTA_BASE_DIR is given, for each species this first looks for
#   FASTA_BASE_DIR/<acc>_results/<acc>_unfiltered.fasta
# and copies it straight in -- no BLAST DB rebuild risk, no internet needed.
# Only if that file is missing for a given species does this fall back to
# `datasets download genome accession <acc> --include protein` (NCBI
# `datasets` CLI, needs outbound internet -- run that fallback on the LOGIN
# node, not via sbatch, since most TACC compute nodes have none).
#
# scer/cneo/calbicans/ncrassa already have proteomes/<short>.faa from
# earlier work and are skipped either way.
#
# Expects species.tsv (short, accession, organism, species_dir) next to it.
# Writes ../proteomes/<short>.faa

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
FASTA_BASE_DIR="${1:-}"
OUT_DIR=../proteomes
mkdir -p "$OUT_DIR"

tail -n +2 species.txt | while IFS=$'\t' read -r short acc organism species_dir; do
    faa="$OUT_DIR/${short}.faa"
    if [[ -s "$faa" ]]; then
        echo "skip: $short ($acc) -- $faa already exists"
        continue
    fi

    if [[ -n "$FASTA_BASE_DIR" ]]; then
        src="$FASTA_BASE_DIR/${acc}_results/${acc}_unfiltered.fasta"
        if [[ -s "$src" ]]; then
            cp "$src" "$faa"
            echo "copied: $short ($acc) <- $src ($(grep -c '^>' "$faa") sequences)"
            continue
        fi
        echo "  ! $src not found for $short, falling back to NCBI download" >&2
    fi

    command -v datasets >/dev/null || {
        echo "ERROR: no unfiltered.fasta found for $short and NCBI 'datasets' CLI" >&2
        echo "       not found on PATH either (check 'module spider datasets')." >&2
        exit 1
    }
    echo "downloading: $short ($acc, $organism) ..."
    tmp=$(mktemp -d)
    (
        cd "$tmp"
        datasets download genome accession "$acc" --include protein
        unzip -q ncbi_dataset.zip
        protein_file=$(find ncbi_dataset/data -name "protein.faa" | head -1)
        if [[ -z "$protein_file" ]]; then
            echo "  ERROR: no protein.faa found for $acc in the downloaded archive" >&2
            exit 1
        fi
        cp "$protein_file" "$OLDPWD/$faa"
    )
    rm -rf "$tmp"
    echo "  wrote $faa ($(grep -c '^>' "$faa") sequences)"
done

echo "done. proteomes present:"
for f in "$OUT_DIR"/*.faa; do
    printf "  %-20s %8d sequences\n" "$(basename "$f" .faa)" "$(grep -c '^>' "$f")"
done
