#!/usr/bin/env Rscript
# Quantify how much of each pool's pre-filtered viral signal the curated panel
# explains: per pool, the fraction of reads mapped on-panel vs. left off-panel
# (unmapped). This consumes only artifacts already produced upstream
# (read_hits.tsv + pool_stats.tsv) and does not depend on R/20. Run from the
# repository root, after R/10_parse_taxonomy.R.
#
#   Rscript R/15_panel_recovery.R
#
# Outputs:
#   results/tables/panel_recovery.tsv        (NEW table; no existing schema changed)
#   results/figures/panel_recovery.png       (stacked on-panel vs off-panel barplot)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

stopifnot(dir.exists("data"), dir.exists("results"))

stats_file <- "results/tables/pool_stats.tsv"
hits_file  <- "results/tables/read_hits.tsv"
if (!file.exists(stats_file)) stop("Missing ", stats_file, " - run scripts/02_seqkit_stats.sh first.")
if (!file.exists(hits_file))  stop("Missing ", hits_file,  " - run R/10_parse_taxonomy.R first.")

# --- total reads per pool from seqkit stats (file basename -> pool_id) -----------
stats <- fread(stats_file)
stats[, pool_id := sub("\\.fasta$", "", basename(file))]
totals <- stats[, .(pool_id, total_reads = num_seqs)]

# --- mapped-read breakdown per pool from read_hits.tsv ---------------------------
# read_hits has one row per mapped read; reads with no panel alignment never
# appear here, so off-panel/unmapped reads are inferred from the pool total.
hits <- fread(hits_file)
mapped <- hits[, .(
  reads_on_panel             = sum(virus != "unassigned_in_panel"),
  reads_unassigned_in_panel  = sum(virus == "unassigned_in_panel")
), by = pool_id]

# Keep every pool present in pool_stats, even pools with zero mapped reads.
rec <- merge(totals, mapped, by = "pool_id", all.x = TRUE)
rec[is.na(reads_on_panel),            reads_on_panel := 0L]
rec[is.na(reads_unassigned_in_panel), reads_unassigned_in_panel := 0L]

# Off-panel = total minus everything that mapped to the panel (any virus, incl.
# the panel-accession-but-unnamed bucket, which is expected ~0).
rec[, reads_offpanel_unmapped := total_reads - (reads_on_panel + reads_unassigned_in_panel)]

rec[, pct_on_panel := round(100 * reads_on_panel / total_reads, 2)]
rec[, pct_offpanel := round(100 * reads_offpanel_unmapped / total_reads, 2)]

# Sanity: the three mutually-exclusive buckets must reconstruct the pool total.
bad <- rec[reads_on_panel + reads_unassigned_in_panel + reads_offpanel_unmapped != total_reads]
if (nrow(bad) > 0L) {
  stop("Per-pool read counts do not sum to total_reads for: ",
       paste(bad$pool_id, collapse = ", "))
}

setcolorder(rec, c("pool_id", "total_reads", "reads_on_panel", "pct_on_panel",
                   "reads_unassigned_in_panel", "reads_offpanel_unmapped", "pct_offpanel"))
setorder(rec, pool_id)

out <- "results/tables/panel_recovery.tsv"
fwrite(rec, out, sep = "\t")
cat(sprintf("Wrote %s (%d pools)\n", out, nrow(rec)))

# --- stacked barplot: on-panel vs off-panel fraction per pool --------------------
long <- melt(
  rec[, .(pool_id,
          on_panel        = reads_on_panel,
          unassigned      = reads_unassigned_in_panel,
          offpanel_unmapped = reads_offpanel_unmapped)],
  id.vars = "pool_id", variable.name = "bucket", value.name = "reads")
long[, bucket := factor(bucket,
  levels = c("on_panel", "unassigned", "offpanel_unmapped"),
  labels = c("on-panel", "unassigned (panel)", "off-panel / unmapped"))]

p <- ggplot(long, aes(x = pool_id, y = reads, fill = bucket)) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_fill_manual(values = c("on-panel"             = "#2c7fb8",
                               "unassigned (panel)"   = "#fdae6b",
                               "off-panel / unmapped" = "#bdbdbd"),
                    name = NULL) +
  labs(title = "Panel recovery: fraction of pre-filtered viral reads explained by the panel",
       x = "pool", y = "fraction of reads") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

fig <- "results/figures/panel_recovery.png"
ggsave(fig, p, width = 10, height = 6, dpi = 150)
cat(sprintf("Wrote %s\n", fig))
