#!/usr/bin/env bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need python
need seqkit
[[ -s "$GENEDIR/download_status.tsv" ]] || die 'Run step 01 first.'
[[ ! -e "$GENEDIR/${GENE}_longest_refseq.fasta" ]] || die 'Extraction output already exists. Move it aside before rerunning.'
python "$ROOT/scripts/pipeline.py" extract "$GENEDIR" "$GENE"
seqkit stats "$GENEDIR/${GENE}_longest_refseq.fasta"
printf 'Review %s/selection.tsv before alignment.\n' "$GENEDIR"
