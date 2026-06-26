# Orchestration for the Zephyr biosurveillance pipeline.
# Run inside the activated conda env (see environment.yml). Targets are ordered.
SHELL := /bin/bash

.PHONY: all pools refs stats taxonomy coverage calls cluster recovery memo clean

all: stats taxonomy coverage calls cluster memo

pools:
	bash scripts/00_fetch_pools.sh

refs:
	bash scripts/01_fetch_refs.sh

stats: pools
	bash scripts/02_seqkit_stats.sh

taxonomy: pools refs
	bash scripts/03_taxonomy_minimap2.sh
	Rscript R/10_parse_taxonomy.R

coverage: pools refs
	bash scripts/04_coverage.sh
	Rscript R/30_coverage_plots.R

# Confidence calls need genome breadth from coverage, so run after both taxonomy and coverage.
calls: taxonomy coverage
	Rscript R/20_confidence_calls.R

cluster: pools
	Rscript R/40_kmer_clustering.R

# Additive: panel-recovery table (depends on read_hits.tsv from `taxonomy`).
recovery: taxonomy
	Rscript R/15_panel_recovery.R

memo:
	Rscript -e 'rmarkdown::render("memo/memo.Rmd")'
	Rscript -e 'rmarkdown::render("memo/memo.Rmd", output_format = "word_document")'

clean:
	rm -rf results/paf/* results/bam/* results/coverage/* \
	       results/tables/* results/figures/* memo/memo.html memo/memo.docx
