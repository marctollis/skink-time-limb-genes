#!/usr/bin/env bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need macse
need python
input="$GENEDIR/${GENE}_longest_refseq.fasta"
[[ -s "$input" ]] || die 'Run step 02 first.'
[[ ! -e "$GENEDIR/alignment" ]] || die 'alignment already exists. Move it aside before rerunning.'
python "$ROOT/scripts/pipeline.py" check "$input" 2
mkdir -p "$GENEDIR/alignment"
macse -prog alignSequences -seq "$input" \
    -out_NT "$GENEDIR/alignment/${GENE}_NT.fasta" \
    -out_AA "$GENEDIR/alignment/${GENE}_AA.fasta" \
    > "$GENEDIR/alignment/macse.log" 2>&1
# Preserve the original MACSE alignment for inspection. Mask whole affected
# codons (including terminal stops), never delete individual characters.
macse -prog exportAlignment -align "$GENEDIR/alignment/${GENE}_NT.fasta" \
    -codonForInternalStop NNN -codonForFinalStop NNN \
    -codonForInternalFS NNN -charForRemainingFS - \
    -out_NT "$GENEDIR/alignment/${GENE}_clean_NT.fasta" \
    > "$GENEDIR/alignment/export.log" 2>&1
python "$ROOT/scripts/pipeline.py" check "$GENEDIR/alignment/${GENE}_clean_NT.fasta" 2 aligned
printf 'Inspect raw NT/AA alignments and logs for frameshifts (!) and stops (*). Cleaned alignment masks these for the training DNA tree.\n'
