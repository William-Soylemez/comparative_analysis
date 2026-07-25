#!/bin/bash
# Download protein FASTA for every species in species.tsv that doesn't
# already have one under proteomes/<short>.faa. Run this on the LOGIN node
# (not via sbatch) -- compute nodes on most TACC systems have no outbound
# internet access, and this just needs the NCBI `datasets` CLI plus a few
# seconds per species.
#
# scer/cneo/calbicans/ncrassa already have proteomes/<short>.faa from
# earlier work; this only fetches the 6 new ones (bdendro, rhizophagus,
# kickxella, agaricus, blastocladiella, mucor).
#
# Usage: bash 01_download_proteomes.sh
# Expects species.tsv (short, accession, organism, species_dir) next to it.
# Writes ../proteomes/<short>.faa

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
OUT_DIR=../proteomes
mkdir -p "$OUT_DIR"

command -v datasets >/dev/null || {
    echo "ERROR: NCBI 'datasets' CLI not found on PATH. On TACC this is usually" >&2
    echo "available via a module (check 'module spider datasets' / 'module spider ncbi-datasets')." >&2
    exit 1
}

tail -n +2 species.tsv | while IFS=$'\t' read -r short acc organism species_dir; do
    faa="$OUT_DIR/${short}.faa"
    if [[ -s "$faa" ]]; then
        echo "skip: $short ($acc) -- $faa already exists"
        continue
    fi

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
