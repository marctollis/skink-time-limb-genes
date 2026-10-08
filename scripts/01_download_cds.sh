#!/usr/bin/env bash
SKINK_LIST_MODE=1
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need python
if [[ -z "${NCBI_EMAIL:-}" ]]; then
    read -r -p 'Contact email for NCBI requests: ' NCBI_EMAIL || die 'No email supplied.'
    export NCBI_EMAIL
fi
python "$ROOT/scripts/ncbi_cds.py" download-list "$ROOT" "${1:-$ROOT/genes.txt}"
