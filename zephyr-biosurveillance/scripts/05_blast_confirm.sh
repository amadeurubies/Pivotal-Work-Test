#!/usr/bin/env bash
# OPTIONAL sanity check: BLAST a small random subsample of reads from one pool against
# NCBI nt (remote). Use only to spot-confirm top calls or check for off-panel viruses.
# Remote BLAST is slow and rate-limited, so keep N small. Requires internet.
#
# Usage: scripts/05_blast_confirm.sh <pool_id> [n_reads] [mode]
#   mode = all       (default) random subsample of the whole pool -> *.blast.tsv
#   mode = offpanel  subsample only reads that did NOT map to the panel (need
#                    results/paf/<pool_id>.paf) -> *.blast_offpanel.tsv. Use this to
#                    identify viral signal the curated panel cannot name (pair with
#                    high pct_offpanel pools from R/15_panel_recovery.R).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pool_id="${1:-}"
n="${2:-25}"
mode="${3:-all}"
[[ -n "$pool_id" ]] || { echo "Usage: $0 <pool_id> [n_reads] [mode: all|offpanel]" >&2; exit 1; }
case "$mode" in
  all|offpanel) ;;
  *) echo "Unknown mode '$mode' (expected: all | offpanel)" >&2; exit 1 ;;
esac

fa="$ROOT/data/raw/${pool_id}.fasta"
[[ -s "$fa" ]] || { echo "Missing $fa" >&2; exit 1; }

sub="$(mktemp --suffix=.fasta)"
cleanup() { rm -f "$sub" "${ids:-}" "${unmapped:-}"; }
trap cleanup EXIT

if [[ "$mode" == "offpanel" ]]; then
  out="$ROOT/results/tables/${pool_id}.blast_offpanel.tsv"
  [[ -s "$out" ]] && { echo "Exists, skipping: $out"; exit 0; }
  paf="$ROOT/results/paf/${pool_id}.paf"
  [[ -s "$paf" ]] || { echo "Missing $paf (run scripts/03_taxonomy_minimap2.sh first)" >&2; exit 1; }
  # Reads that mapped to the panel (PAF col 1); exclude them to keep only off-panel reads.
  ids="$(mktemp)"
  cut -f1 "$paf" | sort -u > "$ids"
  unmapped="$(mktemp --suffix=.fasta)"
  seqkit grep -v -f "$ids" "$fa" > "$unmapped"
  # Reproducible random subsample of n off-panel reads (fixed seed).
  seqkit sample -s 11 -n "$n" "$unmapped" > "$sub"
else
  out="$ROOT/results/tables/${pool_id}.blast.tsv"
  # Reproducible random subsample of n reads (fixed seed).
  seqkit sample -s 11 -n "$n" "$fa" > "$sub"
fi

echo "[blastn -remote] $pool_id ($n reads, mode=$mode) -> $out"
blastn -query "$sub" -db nt -remote \
  -max_target_seqs 1 -evalue 1e-10 \
  -outfmt '6 qseqid sacc stitle pident length evalue bitscore' \
  > "$out"
echo "Wrote $out"
