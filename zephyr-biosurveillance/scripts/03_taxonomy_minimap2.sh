#!/usr/bin/env bash
# Map each pool's reads against the reference panel (ONT preset) and emit one PAF per pool.
# PAF gives, per alignment: residue matches (col 10), alignment block length (col 11),
# and mapping quality (col 12) -> used downstream for identity and confidence.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/raw"
PANEL="$ROOT/data/refs/panel.fasta"
PAFDIR="$ROOT/results/paf"
THREADS="${THREADS:-4}"
mkdir -p "$PAFDIR"

[[ -s "$PANEL" ]] || { echo "Missing panel: $PANEL. Run 01_fetch_refs.sh first." >&2; exit 1; }

shopt -s nullglob
fastas=("$RAW"/*.fasta)
shopt -u nullglob
[[ ${#fastas[@]} -gt 0 ]] || { echo "No FASTAs in $RAW. Run 00_fetch_pools.sh first." >&2; exit 1; }

for fa in "${fastas[@]}"; do
  pool_id="$(basename "$fa" .fasta)"
  out="$PAFDIR/${pool_id}.paf"
  echo "[minimap2] $pool_id"
  # -x map-ont: ONT long-read preset; -c: base-level alignment (accurate match counts);
  # --secondary=no: keep one primary hit per read for clean per-read assignment.
  minimap2 -x map-ont -c --secondary=no -t "$THREADS" "$PANEL" "$fa" > "$out"
done

echo "PAFs in $PAFDIR"
