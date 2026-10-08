#!/usr/bin/env bash
SKINK_LIST_MODE=1
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
need python
need datasets
python "$ROOT/scripts/ncbi_cds.py" download-list "$ROOT" "${1:-$ROOT/genes.txt}"
