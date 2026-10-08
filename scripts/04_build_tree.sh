#!/usr/bin/env bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need iqtree2
need python
input="$GENEDIR/alignment/${GENE}_clean_NT.fasta"
[[ -s "$input" ]] || die 'Run step 03 first.'
threads=${SLURM_CPUS_PER_TASK:-1}
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die 'SLURM_CPUS_PER_TASK must be a positive integer.'
python "$ROOT/scripts/pipeline.py" check "$input" 4 aligned
[[ ! -e "$GENEDIR/iqtree" ]] || die 'iqtree already exists. Move it aside before rerunning.'
mkdir -p "$GENEDIR/iqtree"
iqtree2 -s "$input" -st DNA -m MFP -B 1000 -alrt 1000 \
    -T "$threads" -seed 12345 -pre "$GENEDIR/iqtree/$GENE" \
    > "$GENEDIR/iqtree/console.log" 2>&1
printf 'Tree: %s/iqtree/%s.treefile\nUpdate metadata/gene_tracking.tsv after reviewing the tree.\n' "$GENEDIR" "$GENE"
