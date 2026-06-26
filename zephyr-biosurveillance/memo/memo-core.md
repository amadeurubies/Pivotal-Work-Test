# Pivotal Work Test — Amadeu Rubiés — 26/06/2026

## Time Taken + LLM Usage

4h15m. I didn't know what to expect so had to set up the appropriate environment on my PC +
resubscribe to Cursor, which took about 45m extra at the beginning.

Part 1 is all me with some research. Write-up of the rest is mostly AI for the technical details, etc.

## Part 1: Conceptual Question

- Get a set of real sequences (e.g. from NCBI) and AI-generated sequences from multiple generation
  tools.
- Data matching based on biological characteristics (e.g. GC-content and k-mer frequencies, based on
  preliminary research.)
- If possible, split test data into in-distribution and out-of-distribution sets. OOD set populated
  with more obscure sequences or a non-public database, and those generated from a newly released
  model).
- Build a dumb model as a baseline, e.g. a logistic regression on GC-content and sequence length.
- Evals: False Positive Rate, AUROC, AUPRC, and F-1 score.
- Determine decision thresholds: Depends on the purpose of the evaluation and available resources for
  human review. Assuming this is in a biosecurity context, with a fixed review capacity of N
  sequences/day, set the probability threshold for human review such that (in expectation) N
  sequences are above the threshold per day.

## Data and methodology for Parts 2 and 3

*Full figures and per-pool tables are in the repo: the extended memo (`memo/memo-extended.html`) and
`results/tables/`.*

15 pooled Boston nasal-swab pools, untargeted ONT metatranscriptomics. The public per-pool FASTAs are
reads already classified as vertebrate-infecting viruses. Read counts ≠ prevalence. Reads were mapped
with `minimap2 -x map-ont` to a respiratory-virus RefSeq panel. Per-(pool, virus) calls combine read
support, median identity, and genome breadth into high/medium/low tiers (high = ≥50 reads, identity
≥0.90, breadth ≥25%). Track B embeds reads as tetranucleotide (k=4) frequencies → UMAP → clustering,
scored against the minimap2 labels by adjusted Rand index (ARI).

## Part 2: Taxonomic Classification

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

- **Seasonal coronaviruses are highest prevalence.** HCoV-229E dominates five pools (≈78k–100k reads,
  75–100% breadth, ≈8,000–12,000× depth); HCoV-NL63 dominates pool11; HCoV-HKU1 gives a near-complete
  genome in pool03. At ONT-typical identity (≈0.97–0.99) these are unambiguous, near-complete genomes
  and most high-confidence calls.
- **Rhinovirus is ubiquitous.** Its "low" tier is a reference artefact rather than weak evidence for
  its prevalence. RV-A/B/C are often the most abundant taxon (e.g. pool10: 96,854 reads covering
  99.7% of the ~7 kb genome at high depth), yet they fall in the "low" tier only because they map at
  ~0.69–0.78 identity. That low identity reflects rhinovirus's high genetic diversity against our
  single RefSeq genome. Reads of the circulating strain match a divergent reference poorly even though
  they are clearly rhinovirus. So "low" means "no close reference in the panel," not "probably absent."
- **Other viruses appear sporadically, and how much of a pool the panel explains tracks the dominant
  taxon.** Besides the seasonal coronaviruses and rhinovirus we see one solid SARS-CoV-2 call
  (pool08), strong hMPV in pool15, and low-to-medium influenza A/B. Where members of the same family
  co-occur in a pool they are flagged as possible cross-mappings. Coverage of the dominant viruses is
  near-complete, and coronaviruses show a 5′/3′-skewed depth profile — apparently since ONT
  metatranscriptomics captures their nested subgenomic mRNAs. The reference panel accounts for >99% of
  reads in coronavirus-dominated pools but only ~34–36% in pools where divergent rhinovirus or
  off-panel viruses dominate (pools 13 and 15), so the unexplained reads are exactly where divergent
  or unpanelled viruses are expected.

## Part 3: Track B Embedding-based Clustering

*Figure 1 — k-mer UMAP (`results/figures/kmer_umap_taxonomy.png`).*

K-mer (tetranucleotide) composition reproduces the broad taxonomy but not species-level detail:
ARI ~0.42 (hierarchical), 0.40 (k-means, k=8), 0.05 (DBSCAN). The modest score is mostly
fragmentation rather than cross-taxon mixing — coronaviruses and rhinoviruses occupy separate UMAP
regions and never share a cluster, but a single taxon is spread across several clusters (HCoV-229E
across four) and the three rhinovirus species (A/B/C) collapse together. So the clusters corroborate
the alignment-based calls at the genus level. The unassigned reads never form a separate group and
never join the coronavirus clusters. They fall entirely within the rhinovirus clusters, consistent
with the Part 2 reading that they are divergent rhinovirus/picornavirus below the identity gate rather
than a novel taxon.

## Limitations

100k-read/pool subsample, for time reasons. Reads are pre-filtered and pooled. Off-panel reads are
unassigned (a remote-BLAST spot-check of the high-off-panel pools is implemented but didn't finish in
time). ONT error (~5–15%) lowers identity, so identity alone never calls a virus and close relatives
that may cross-map are flagged.

## Reproducibility

GitHub. Single conda environment (`environment.yml`); fixed run order in the repo README/`Makefile`
(confidence calls run after coverage); seeds set for UMAP/k-means. The committed results use a
100k-read/pool subsample (`MAX_READS=100000 THREADS=12`, seed 11); the full-data run uses the same
commands without `MAX_READS`.
