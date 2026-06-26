# Agent instructions

## Project overview

- **Purpose:** R-first metagenomic biosurveillance of SecureBio Zephyr pooled nasal-swab viral-read
  FASTAs (ONT long reads). Confirms which viruses are present and with what confidence, assesses
  genome coverage (Part 2), and clusters reads by k-mer composition vs taxonomy (Part 3, Track B).
- **Main stack:** R (`data.table`, `Biostrings`, `uwot`, `dbscan`, `ggplot2`, `rmarkdown`) plus
  Unix bioinformatics tools (`minimap2`, `samtools`, `bedtools`, `seqkit`, BLAST+), all pinned in a
  single conda/mamba environment (`environment.yml`).
- **Important directories** (under `zephyr-biosurveillance/`): `scripts/` (shell fetch/align/coverage
  steps), `R/` (parsing, confidence calls, coverage plots, k-mer clustering), `data/` (inputs;
  `raw/` and large refs are gitignored), `results/` (tables + figures), `memo/` (write-up).
- **Runtime assumptions:** runs in **WSL2 / Linux** with the `zephyr` conda env active. The pipeline
  does **not** run on native Windows. Tools live inside the env, not on the base `PATH`.

## Commands

Run from `zephyr-biosurveillance/` with the env active (`conda activate zephyr`):

- **Install:** `mamba env create -f environment.yml` (or `conda env create -f environment.yml`)
- **Run full pipeline:** `make all` (see README "Run order" for the step-by-step equivalent)
- **Render memo:** `Rscript -e 'rmarkdown::render("memo/memo-extended.Rmd")'`
  (add `output_format = "word_document"` for the `.docx`)
- **Test / Typecheck / Lint / Build:** none — this is an analysis pipeline, not an application.
  Validate by running the relevant pipeline step and checking the outputs in `results/`.

## Working rules

- Inspect relevant files before editing.
- Keep changes focused.
- Follow existing project conventions (R-first, shell tools in one conda env, fixed run order).
- There is no automated test suite; verify behaviour changes by re-running the affected step and
  inspecting `results/` tables/figures.
- Do not edit generated, vendored, secret, or production config files unless explicitly asked.
- Do not commit raw FASTAs, large intermediates, or rendered memos (see `.gitignore`).
- Do not introduce dependencies without justification.
- Summarise changes and validation after each task.
