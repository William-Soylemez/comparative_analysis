#!/bin/bash
# Step 1/3: reciprocal BLASTP for a pair of species.
#
# For each species, build a protein BLAST DB from proteomes/<short>.faa (if not
# already built), then run blastp in BOTH directions. The bitscore column of
# these hits is the raw signal that 2_similarity.py turns into the IsoRank
# sequence-similarity input.
#
# Usage: bash 1_blast.sh <a> <b>
#   e.g. bash 1_blast.sh scer calbicans
# Expects: proteomes/<a>.faa, proteomes/<b>.faa, and blastp + makeblastdb on PATH.
# Writes:  blastdb/<short>.*, blast_out/<a>_vs_<b>.tsv, blast_out/<b>_vs_<a>.tsv
#          (tabular outfmt 6: qseqid sseqid bitscore evalue pident)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

A="${1:?usage: bash 1_blast.sh <a> <b>}"
B="${2:?usage: bash 1_blast.sh <a> <b>}"
EVALUE="${EVALUE:-1e-5}"
FMT="6 qseqid sseqid bitscore evalue pident"

for tool in makeblastdb blastp; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not found on PATH" >&2; exit 1; }
done

mkdir -p blastdb blast_out

build_db() {
    local s="$1"
    [[ -s "blastdb/${s}.pin" ]] && { echo "skip db: $s"; return 0; }
    [[ -s "proteomes/${s}.faa" ]] || { echo "ERROR: missing proteomes/${s}.faa" >&2; exit 1; }
    echo "building blastdb: $s"
    makeblastdb -in "proteomes/${s}.faa" -dbtype prot -out "blastdb/${s}" >/dev/null
}

run_blast() {
    local q="$1" s="$2" out="blast_out/${1}_vs_${2}.tsv"
    if [[ -s "$out" ]]; then echo "skip blast: ${q}_vs_${s}"; return 0; fi
    blastp -query "proteomes/${q}.faa" -db "blastdb/${s}" \
           -evalue "$EVALUE" -outfmt "$FMT" -out "$out.tmp"
    mv "$out.tmp" "$out"
    echo "done blast: ${q}_vs_${s} ($(wc -l < "$out") hits)"
}

build_db "$A"
build_db "$B"
run_blast "$A" "$B"
run_blast "$B" "$A"
