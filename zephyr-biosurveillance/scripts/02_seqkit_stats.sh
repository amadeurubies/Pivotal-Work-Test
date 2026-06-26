#!/usr/bin/env bash
# Per-pool read statistics (counts, length distribution, N50) for all downloaded pools.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"
OUT="$ROOT/results/tables/pool_stats.tsv"
mkdir -p "$ROOT/results/tables"

shopt -s nullglob
fastas=("$RAW"/*.fasta)
shopt -u nullglob
[[ ${#fastas[@]} -gt 0 ]] || { echo "No FASTAs in $RAW. Run 00_fetch_pools.sh first." >&2; exit 1; }

# -a: all stats (incl. N50, min/max len); -T: tab-separated output.
seqkit stats -a -T "${fastas[@]}" > "$OUT"
echo "Wrote $OUT"
