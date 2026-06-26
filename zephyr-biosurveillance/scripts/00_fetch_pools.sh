#!/usr/bin/env bash
# Download per-pool viral-read FASTAs listed in data/metadata/pools.tsv.
# Fill the `url` column in pools.tsv first (copy the FASTA link from the
# Zephyr sample log: https://data.securebio.org/zephyr/#respiratory-viral-reads).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$ROOT/data/metadata/pools.tsv"
OUTDIR="$ROOT/data/raw"
# Optional speed lever: cap reads per pool. Empty/unset = no cap = full data (default).
# Capping happens here so every downstream step inherits the smaller set automatically.
MAX_READS="${MAX_READS:-}"
SAMPLE_SEED=11   # fixed seed: downsampling must be reproducible (repo requires seeds)
mkdir -p "$OUTDIR"

[[ -f "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 1; }

n_ok=0
n_skip=0
# Skip header; read tab-separated columns.
while IFS=$'\t' read -r pool_id date site pool_size n_reads url; do
  # Strip any stray carriage returns (in case the file was edited on Windows).
  pool_id="${pool_id%$'\r'}"; url="${url%$'\r'}"
  [[ -z "${pool_id:-}" || "$pool_id" == "pool_id" ]] && continue

  if [[ -z "${url:-}" || "$url" == TODO_* ]]; then
    echo "[skip] $pool_id: url not set in manifest" >&2
    n_skip=$((n_skip + 1))
    continue
  fi

  dest="$OUTDIR/${pool_id}.fasta"
  if [[ -s "$dest" ]]; then
    echo "[have] $dest"
    n_ok=$((n_ok + 1))
    continue
  fi

  echo "[get ] $pool_id <- $url"
  tmp="$(mktemp)"
  if ! curl -fSL "$url" -o "$tmp"; then
    echo "[FAIL] $pool_id: download error" >&2
    rm -f "$tmp"
    continue
  fi
  # Decompress if the payload is gzip, otherwise use as-is.
  raw="$(mktemp)"
  if gzip -t "$tmp" 2>/dev/null; then
    gzip -dc "$tmp" > "$raw"
  else
    cp "$tmp" "$raw"
  fi
  rm -f "$tmp"

  # Apply the optional read cap. With MAX_READS unset the file is moved verbatim
  # (byte-for-byte the decompressed payload). When set, pools above the cap are
  # downsampled to MAX_READS reads with a fixed seed; pools at/below it are kept as-is.
  if [[ -n "$MAX_READS" ]] && (( $(seqkit stats -T "$raw" | awk 'NR==2 {print $4}') > MAX_READS )); then
    echo "[cap ] $pool_id -> $MAX_READS reads (seed $SAMPLE_SEED)"
    seqkit sample -s "$SAMPLE_SEED" -n "$MAX_READS" -2 "$raw" > "$dest"
    rm -f "$raw"
  else
    mv "$raw" "$dest"
  fi
  n_ok=$((n_ok + 1))
done < <(tail -n +2 "$MANIFEST")

echo "Pools ready: $n_ok | skipped (no url): $n_skip"
echo "FASTAs in $OUTDIR"
