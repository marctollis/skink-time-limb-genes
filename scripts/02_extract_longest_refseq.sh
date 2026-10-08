#!/usr/bin/env bash
# Select longest eligible RefSeq CDS per species from the clade-wide download.
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need python
need seqkit
python "$ROOT/scripts/ncbi_cds.py" extract "$GENEDIR" "$GENE"
seqkit stats "$GENEDIR/${GENE}_longest_cds.fasta"
printf 'Review candidate_cds.tsv and selection.tsv in %s before alignment.\n' "$GENEDIR"
