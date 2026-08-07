#!/bin/bash
# Build a BLAST protein DB under blastdb/<short>.* for every species with a
# proteomes/<short>.faa that doesn't already have one. Cheap (seconds per
# species) -- fine to run on the login node right after 01_download_proteomes.sh,
# or as a quick non-parallel step at the top of the sbatch blast job.
#
# Usage: bash 02_build_blastdbs.sh
# Expects ../species.txt, ../proteomes/<short>.faa
# Writes ../blastdb/<short>.{phr,pin,psq,...}

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p ../blastdb

command -v makeblastdb >/dev/null || {
    echo "ERROR: makeblastdb not found on PATH. Vista has no 'blast' module" \
         "(confirmed via 'module spider blast' -> not found) -- if running this" \
         "standalone (not via 00_run_all.sbatch/03_run_all_pairs_blast.sbatch," \
         "which already set this), first run:" \
         "  export PATH=\"/work/11301/wsoylemez/vista/ncbi-blast-current/bin:\$PATH\"" >&2
    exit 1
}

tail -n +2 ../species.txt | while IFS=$'\t' read -r short acc organism; do
    if [[ -s "../blastdb/${short}.pin" ]]; then
        echo "skip: $short -- blastdb/${short}.pin already exists"
        continue
    fi
    faa="../proteomes/${short}.faa"
    if [[ ! -s "$faa" ]]; then
        echo "  ! skipping $short: no proteomes/${short}.faa (run 01_download_proteomes.sh first)" >&2
        continue
    fi
    echo "building blastdb: $short"
    makeblastdb -in "$faa" -dbtype prot -out "../blastdb/${short}"
done
