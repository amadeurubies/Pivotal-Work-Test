#!/usr/bin/env bash
# Align each pool to the panel, then compute genome coverage against references.
# Outputs per pool:
#   results/coverage/<pool>.coverage.txt  - samtools coverage (breadth % + mean depth per reference)
#   results/coverage/<pool>.bedgraph      - bedtools genomecov (depth across genome positions, for plots)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"
PANEL="$ROOT/data/refs/panel.fasta"
BAMDIR="$ROOT/results/bam"
COVDIR="$ROOT/results/coverage"
THREADS="${THREADS:-4}"
mkdir -p "$BAMDIR" "$COVDIR"

[[ -s "$PANEL" ]] || { echo "Missing panel: $PANEL. Run 01_fetch_refs.sh first." >&2; exit 1; }

shopt -s nullglob
fastas=("$RAW"/*.fasta)
shopt -u nullglob
[[ ${#fastas[@]} -gt 0 ]] || { echo "No FASTAs in $RAW. Run 00_fetch_pools.sh first." >&2; exit 1; }

for fa in "${fastas[@]}"; do
  pool_id="$(basename "$fa" .fasta)"
  bam="$BAMDIR/${pool_id}.bam"
  echo "[align] $pool_id"
  # Exclude unmapped (0x4), secondary (0x100) and supplementary (0x800) for honest depth.
  minimap2 -ax map-ont -t "$THREADS" "$PANEL" "$fa" \
    | samtools view -b -F 0x904 - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"

  samtools coverage "$bam" > "$COVDIR/${pool_id}.coverage.txt"
  bedtools genomecov -ibam "$bam" -bga > "$COVDIR/${pool_id}.bedgraph"
done

echo "Coverage outputs in $COVDIR"
