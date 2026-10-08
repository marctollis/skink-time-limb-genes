#!/usr/bin/env bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need datasets
need unzip
need python
[[ -s "$ROOT/species.txt" ]] || die 'species.txt is missing or empty.'
[[ ! -e "$GENEDIR/downloads" ]] || die 'downloads already exists. Use a new gene, or move the old gene directory aside before starting again.'
python "$ROOT/scripts/pipeline.py" species "$ROOT/species.txt" > "$GENEDIR/species.tsv"
mkdir -p "$GENEDIR/downloads"
printf 'species\tlabel\tstatus\n' > "$GENEDIR/download_status.tsv"
failed=0
while IFS=$'\t' read -r species label; do
    dest="$GENEDIR/downloads/$label"
    mkdir -p "$dest"
    printf 'Downloading %s / %s\n' "$GENE" "$species"
    if datasets download gene symbol "$GENE" --taxon "$species" --include cds --filename "$dest/package.zip" > "$dest/download.log" 2>&1; then
        if unzip -q "$dest/package.zip" -d "$dest/package" >> "$dest/download.log" 2>&1; then
            printf '%s\t%s\tdownloaded\n' "$species" "$label" >> "$GENEDIR/download_status.tsv"
        else
            printf '%s\t%s\tunzip_failed\n' "$species" "$label" >> "$GENEDIR/download_status.tsv"
            failed=1
        fi
    else
        printf '%s\t%s\tdownload_failed\n' "$species" "$label" >> "$GENEDIR/download_status.tsv"
        failed=1
    fi
done < "$GENEDIR/species.tsv"
printf 'Review %s/download_status.tsv and per-species download.log files.\n' "$GENEDIR"
(( failed == 0 )) || die 'One or more requests failed (including potentially unavailable genes). Successful packages were retained; review failures before step 02.'
