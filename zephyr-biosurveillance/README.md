# Zephyr metagenomic biosurveillance

Lightweight, R-first analysis of SecureBio **Zephyr** pooled nasal-swab viral-read FASTAs
(ONT long reads, pre-filtered to vertebrate-infecting viruses). It answers, across >=10 pools:
which viruses are present and with what confidence, what genome coverage we see for common
respiratory viruses (Part 2), and how reads cluster by k-mer composition vs taxonomy (Part 3,
Track B).

> Designed to run in **WSL2 (Ubuntu)** with a single **conda/mamba** environment. The core
> tools (`minimap2`, `samtools`, `bedtools`, `seqkit`, BLAST+) are Unix-native.

## Requirements

```bash
# In WSL2 / Linux, with conda or mamba installed:
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate zephyr
```

> **First-time setup notes (WSL2):**
> - If you don't already have `conda`/`mamba`, install a base distribution first (e.g.
>   [Miniforge](https://github.com/conda-forge/miniforge)), then run the commands above.
> - **All run-order commands below assume the `zephyr` env is active.** The pipeline tools
>   (`seqkit`, `minimap2`, `samtools`, `bedtools`, BLAST+, `Rscript`) live *inside* the env and
>   are **not** on the base `PATH` — running a script before `conda activate zephyr` will fail
>   with "command not found". Sanity-check with `command -v seqkit minimap2 Rscript`.
> - This pipeline does **not** run on native Windows; run everything from a WSL2 shell.

## Data

The pool manifest [`data/metadata/pools.tsv`](data/metadata/pools.tsv) is pre-filled with the 15
largest pools (date, site, pool size, read count, FASTA URL) from the Zephyr sample log at
`https://data.securebio.org/zephyr/#respiratory-viral-reads`. The `url` column is already
populated; to swap in different pools, copy the link from the relevant read-count cell. The
reference panel is defined in [`data/refs/refs.tsv`](data/refs/refs.tsv) (RefSeq accessions for
common respiratory viruses); edit it if an accession does not resolve.

## Run order

```bash
# 1. Fetch inputs
bash scripts/00_fetch_pools.sh        # pool FASTAs  -> data/raw/
bash scripts/01_fetch_refs.sh         # panel.fasta  -> data/refs/

# 2. Read stats + alignment-based taxonomy
bash scripts/02_seqkit_stats.sh       # -> results/tables/pool_stats.tsv
bash scripts/03_taxonomy_minimap2.sh  # -> results/paf/*.paf
Rscript R/10_parse_taxonomy.R         # -> results/tables/read_hits.tsv
Rscript R/15_panel_recovery.R         # -> panel_recovery.tsv (on- vs off-panel read fractions)

# 3. Coverage, then confidence calls
#    R/20 needs genome breadth, so coverage (04 + R/30) MUST run before it.
bash scripts/04_coverage.sh           # -> results/coverage/*
Rscript R/30_coverage_plots.R         # -> coverage_summary.tsv + plots
Rscript R/20_confidence_calls.R       # -> taxonomy_calls.tsv + heatmap (uses coverage breadth)

# 4. Track B: k-mer embedding + clustering
Rscript R/40_kmer_clustering.R        # -> cluster_* tables + UMAP plots

# 5. Memo
#   - memo/memo-core.md     : concise core memo (Parts 2-3) -> paste into the Google Doc.
#   - memo/memo-extended.Rmd: extended memo (same narrative with all figures + full tables).
Rscript -e 'rmarkdown::render("memo/memo-extended.Rmd")'                                  # -> memo/memo-extended.html
Rscript -e 'rmarkdown::render("memo/memo-extended.Rmd", output_format = "word_document")'  # -> memo/memo-extended.docx

# Optional: remote BLAST spot-check of one pool
bash scripts/05_blast_confirm.sh pool01 25
# Optional: BLAST only the off-panel (unmapped) reads of a high-pct_offpanel pool
bash scripts/05_blast_confirm.sh pool01 25 offpanel   # -> *.blast_offpanel.tsv
```

Or run everything with `make all` (see [`Makefile`](Makefile)).

### Speed levers (optional)

The default flow above uses the full read set and 4 threads, and takes ~1.5-3h on a 6c/12t
machine. Two optional levers can bring a run under an hour without changing default behaviour
or output formats:

- **Threads** — `scripts/03` and `scripts/04` honour a `THREADS` env var (default 4, kept safe
  for low-core machines). To use all logical cores, prefix commands, e.g.
  `THREADS=12 bash scripts/03_taxonomy_minimap2.sh`. `make THREADS=12 all` works too.
- **Subsampling** — `scripts/00` honours an optional `MAX_READS` env var (unset = no cap = full
  data, the default). When set, any pool larger than the cap is downsampled to `MAX_READS` reads
  with a fixed seed (`seqkit sample -s 11 -n "$MAX_READS" -2`); smaller pools are left untouched.
  Because the capped FASTAs land in `data/raw/`, every downstream step inherits the smaller set
  automatically.

**Fast run** (taxonomy confirmation in well under an hour):

```bash
MAX_READS=100000 THREADS=12 bash scripts/00_fetch_pools.sh   # capped FASTAs -> data/raw/
# ...then run scripts/01..05 and the R steps exactly as in "Run order" above,
#    prefixing 03/04 with THREADS=12 (e.g. THREADS=12 make all).
```

Tradeoff: capping lowers depth and breadth, so the fast path is fine for **taxonomy
confirmation** but not for near-complete genomes. The full-data flow above remains the primary
documented path. Downsampling is applied once, at download time, and `scripts/00` skips pools
whose FASTA already exists — so to re-cap at a different `MAX_READS`, clear `data/raw/` first
(`rm -f data/raw/*.fasta`) and re-run.

### How the committed results were generated (subsampled fast path)

The tables and figures currently in `results/` were produced with the **subsampled fast path**:
each of the 15 pools was capped to **100,000 reads** (pool15 was already smaller at 76,508), via
`MAX_READS=100000 THREADS=12` and seed 11, then run through taxonomy → coverage → confidence →
clustering. All read counts, breadth, and depth therefore reflect **subsampled input** and are
intended for taxonomy confirmation and relative comparison, not absolute genome completeness. To
reproduce full-data results, clear `data/raw/` and re-run the pipeline without `MAX_READS`.

## Outputs

| File | Description |
| --- | --- |
| `results/tables/pool_stats.tsv` | Per-pool read counts, lengths, N50 |
| `results/tables/read_hits.tsv` | Best minimap2 hit per read, with taxonomy |
| `results/tables/panel_recovery.tsv` | Per-pool on-panel vs off-panel (unmapped) read fractions |
| `results/tables/taxonomy_calls.tsv` | Per-(pool, virus) calls with confidence tiers |
| `results/tables/coverage_summary.tsv` | Breadth (%) and mean depth per reference |
| `results/tables/cluster_vs_taxonomy.tsv`, `cluster_metrics.tsv` | Cluster/taxonomy crosstab + ARI |
| `results/figures/*.png` | Taxonomy heatmap, coverage plots, k-mer UMAPs |
| `memo/memo-core.md` | Concise core memo (Parts 2-3), copy-paste into the Google Doc |
| `memo/memo-extended.html`, `memo/memo-extended.docx` | Extended memo: same narrative with all figures + full tables |

## Confidence criteria (heuristic, tunable)

Per (pool, virus), from read support, median identity, and genome breadth:

- **High:** >=50 reads, median identity >=0.90, breadth >=25%.
- **Medium:** >=5 reads, identity >=0.80, breadth >=5%.
- **Low:** otherwise.

Thresholds are relaxed for ONT error and are defaults, not ground truth. Co-called members of the
same viral family are flagged (`ambiguity_flag`) for possible cross-mapping.

## Notes / caveats

- Reads are pre-filtered and pooled: confirmation, not discovery; abundance != prevalence.
- Anything outside the reference panel is `unassigned`; use `05_blast_confirm.sh` to spot-check.
- Per-read k-mer clusters are noisy and may reflect read length/composition, not only taxonomy.
- The `n_reads` column in `pools.tsv` is transcribed from the sample log for reference only. When a
  date/site has several pools, the URL's `P1/P2/P3` suffix may not match read-count order, so the
  listed `n_reads` could correspond to a sibling pool. This does not affect results: `seqkit`
  recomputes true per-pool counts (`pool_stats.tsv`) from the downloaded FASTAs.
- For segmented viruses (influenza), per-virus `breadth_pct`/`mean_depth` in `taxonomy_calls.tsv`
  are the best single-segment values (max across segments), so they reflect the best-covered
  segment rather than the whole genome. Per-segment detail is in `coverage_summary.tsv`.
- Confidence calls (`R/20`) require genome breadth, so they must run **after** coverage
  (`04` + `R/30`); the run order above and `make all` are sequenced accordingly (the `calls`
  target). Running `R/20` with no coverage present yields breadth = 0 and collapses every call to
  `low`.
- `R/10` reads the 12 fixed PAF columns explicitly (`cut -f1-12`) so it does not silently drop
  alignments when minimap2 emits a variable number of optional tag columns.

**Submission:** paste the concise core memo `memo/memo-core.md` into the Google Doc and format it
there alongside Part 1 (the conceptual question) and the time taken. Insert the two essential figures
it references (taxonomy heatmap, k-mer UMAP) from `results/figures/`. For the fuller version with
every figure and the full per-pool tables, render the extended memo (`memo/memo-extended.Rmd` ->
`output_format = "word_document"`), upload `memo/memo-extended.docx` to Drive, and "Open with Google
Docs" (browser copy-paste tends to drop the figures).
