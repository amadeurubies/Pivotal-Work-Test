#!/usr/bin/env bash
# Fetch the curated respiratory-virus reference panel listed in data/refs/refs.tsv
# (one RefSeq accession per row) via NCBI efetch, concatenate into panel.fasta,
# and record reference lengths. Accessions that fail to fetch are reported and skipped.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFDIR="$ROOT/data/refs"
REFS_TSV="$REFDIR/refs.tsv"
PANEL="$REFDIR/panel.fasta"
mkdir -p "$REFDIR/seqs"

[[ -f "$REFS_TSV" ]] || { echo "Missing $REFS_TSV" >&2; exit 1; }
command -v efetch >/dev/null || { echo "efetch not found (install entrez-direct)" >&2; exit 1; }

: > "$PANEL"
n_ok=0
n_fail=0
while IFS=$'\t' read -r accession virus segment family; do
  accession="${accession%$'\r'}"
  [[ -z "${accession:-}" || "$accession" == "accession" ]] && continue

  seq="$REFDIR/seqs/${accession}.fasta"
  if [[ ! -s "$seq" ]]; then
    echo "[efetch] $accession ($virus ${segment})"
    if ! efetch -db nuccore -id "$accession" -format fasta > "$seq" 2>/dev/null || [[ ! -s "$seq" ]]; then
      echo "[FAIL] $accession could not be fetched" >&2
      rm -f "$seq"
      n_fail=$((n_fail + 1))
      continue
    fi
  fi
  cat "$seq" >> "$PANEL"
  n_ok=$((n_ok + 1))
done < <(tail -n +2 "$REFS_TSV")

echo "References fetched: $n_ok | failed: $n_fail"

# Record reference sequence lengths (accession <tab> length) for downstream R.
if command -v samtools >/dev/null; then
  samtools faidx "$PANEL"
  cut -f1,2 "${PANEL}.fai" > "$REFDIR/ref_lengths.tsv"
elif command -v seqkit >/dev/null; then
  seqkit fx2tab -nl "$PANEL" | awk -F'\t' '{split($1,a," "); print a[1]"\t"$2}' > "$REFDIR/ref_lengths.tsv"
fi

echo "Panel written to $PANEL"
