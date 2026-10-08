#!/usr/bin/env bash
# Shared by the four numbered entry points; paths work from any directory.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing $1; activate the conda environment."; }
if [[ "${SKINK_LIST_MODE:-0}" != 1 ]]; then
    read -r -p 'Gene symbol (e.g. SHH): ' GENE || die 'No gene supplied.'
    [[ "$GENE" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die 'Use a gene symbol containing letters, numbers, dots, underscores or hyphens.'
    GENE=$(printf '%s' "$GENE" | tr '[:lower:]' '[:upper:]')
    GENEDIR="$ROOT/genes/$GENE"
    mkdir -p "$GENEDIR"
fi
