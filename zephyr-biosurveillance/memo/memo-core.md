*Full figures and per-pool tables are in the repo: the extended memo (`memo/memo-extended.html`,
from `memo/memo-extended.Rmd`) and `results/tables/`.*

**Data and methods.** 15 pooled Boston nasal-swab pools, untargeted ONT metatranscriptomics; the
public per-pool FASTAs are reads already classified as vertebrate-infecting viruses. So this is
taxonomic **confirmation and coverage**, not discovery, and read counts ≠ prevalence. Reads were
mapped with `minimap2 -x map-ont` to a curated respiratory-virus RefSeq panel; per-(pool, virus)
calls combine read support, median identity, and genome breadth into high/medium/low tiers
(high = ≥50 reads, identity ≥0.90, breadth ≥25%). Track B embeds reads as tetranucleotide (k=4)
frequencies → UMAP → clustering, scored against the minimap2 labels by adjusted Rand index (ARI).

## Part 2: Taxonomic Classification

*Figure 1 — taxonomy heatmap (`results/figures/taxonomy_heatmap.png`): reads per pool × virus,
tile border = confidence.*

Top high-confidence calls (13 in total; full table: `results/tables/taxonomy_calls.tsv`):

| pool | virus | reads | identity | breadth % | depth × |
| --- | --- | ---: | ---: | ---: | ---: |
| pool11 | HCoV-NL63 | 99,980 | 0.976 | 88.4 | 9,445 |
| pool06 | HCoV-229E | 99,849 | 0.971 | 99.9 | 11,091 |
| pool02 | HCoV-229E | 99,760 | 0.972 | 99.5 | 12,325 |
| pool09 | HCoV-229E | 97,338 | 0.973 | 96.3 | 10,989 |
| pool08 | HCoV-229E | 90,540 | 0.972 | 100.0 | 8,774 |

*+8 more high-confidence calls (HCoV-229E pool04; HCoV-NL63 pools 03/07/09/12; HCoV-HKU1 pool03;
hMPV pool15; SARS-CoV-2 pool08).*

- **Seasonal coronaviruses are the clearest signal.** HCoV-229E dominates five pools (≈78k–100k
  reads, 75–100% breadth, ≈8,000–12,000× depth); HCoV-NL63 dominates pool11; HCoV-HKU1 gives a
  near-complete genome in pool03. At ONT-typical identity (≈0.97–0.99) these are unambiguous,
  near-complete genomes and most high-confidence calls.
- **Rhinovirus is ubiquitous but scores "low" from reference divergence, not absence.** RV-A/B/C
  are often the most abundant taxon (pool10 RV-A 96,854 reads, 99.7% breadth) yet map at only
  ≈0.69–0.78 identity to the single divergent RefSeq while covering the whole ≈7 kb genome deeply —
  here breadth/depth, not identity, are the evidence.
- **Sporadic co-detections and coverage.** SARS-CoV-2 (one solid call, pool08), hMPV (high in
  pool15), and influenza A/B (low–medium) appear at lower support; co-occurring family members are
  flagged as possible cross-mapping. Coverage is near-complete for dominant viruses (coronaviruses
  show a 5′/3′ profile from subgenomic mRNAs), and the panel explains >99% of reads in CoV-dominated
  pools but only 34–36% where divergent rhinovirus/off-panel viruses prevail (pools 13/15).

## Part 3: Track B Embedding-based Clustering

*Figure 2 — k-mer UMAP (`results/figures/kmer_umap_taxonomy.png`).*

K-mer composition recovers the broad taxonomy but over-splits each taxon: ARI 0.42 (hierarchical),
0.40 (k-means, k=8), 0.05 (DBSCAN). Coronaviruses and picornaviruses sit in distinct regions and
clusters are essentially single-taxon pure, so the modest ARI is **fragmentation, not mixing** —
clustering is confirmatory, not a replacement for alignment. Notably, the unassigned reads do not
form their own group but co-cluster with rhinovirus, supporting the Part 2 reading that they are
divergent rhinovirus below the identity gate rather than a novel taxon.

## Limitations

100k-read/pool subsample (relative patterns, not absolute genome completeness); reads are
pre-filtered and pooled (confirmation, not prevalence); off-panel reads are unassigned (a
remote-BLAST spot-check of the high-off-panel pools is implemented but did not finish within the
time budget); ONT error (≈5–15%) lowers identity, so identity alone never calls a virus and close
relatives that may cross-map are flagged.

## Reproducibility

Single conda environment (`environment.yml`); fixed run order in the repo README/`Makefile`
(confidence calls run after coverage); seeds set for UMAP/k-means. The committed results use a
100k-read/pool subsample (`MAX_READS=100000 THREADS=12`, seed 11); the full-data run uses the same
commands without `MAX_READS`. Confidence thresholds are tunable defaults, not ground truth.
