#!/usr/bin/env Rscript
# Summarise genome coverage for the respiratory panel and plot breadth/depth.
# Run from the repository root, after scripts/04_coverage.sh.
#
#   Rscript R/30_coverage_plots.R
#
# Outputs:
#   results/tables/coverage_summary.tsv
#   results/figures/coverage_breadth.png
#   results/figures/coverage_profiles.png   (depth across genome for top hits)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

stopifnot(dir.exists("data"), dir.exists("results"))

refs <- fread("data/refs/refs.tsv")

cov_files <- list.files("results/coverage", pattern = "\\.coverage\\.txt$", full.names = TRUE)
if (length(cov_files) == 0L) stop("No coverage files - run scripts/04_coverage.sh first.")

read_cov <- function(f) {
  dt <- fread(f)
  setnames(dt, 1, "rname")
  dt[, pool_id := sub("\\.coverage\\.txt$", "", basename(f))]
  dt
}
cov <- rbindlist(lapply(cov_files, read_cov), use.names = TRUE, fill = TRUE)
cov <- merge(cov, refs, by.x = "rname", by.y = "accession", all.x = TRUE)

# Per (pool, virus, segment) summary; keep references with at least one read.
summary <- cov[numreads > 0 & !is.na(virus),
  .(pool_id, virus, segment, accession = rname,
    ref_len = endpos, n_reads = numreads,
    breadth_pct = round(coverage, 2), mean_depth = round(meandepth, 2))]
setorder(summary, pool_id, -n_reads)

out <- "results/tables/coverage_summary.tsv"
fwrite(summary, out, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", out, nrow(summary)))

# --- breadth tile plot (virus x pool) -------------------------------------------
breadth_by_virus <- cov[!is.na(virus),
  .(breadth_pct = max(coverage, na.rm = TRUE)), by = .(pool_id, virus)]
pb <- ggplot(breadth_by_virus, aes(pool_id, virus, fill = breadth_pct)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "breadth (%)", limits = c(0, 100)) +
  labs(title = "Genome breadth of coverage (>=1x) by pool", x = "pool", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave("results/figures/coverage_breadth.png", pb, width = 10, height = 6, dpi = 150)
cat("Wrote results/figures/coverage_breadth.png\n")

# --- depth profiles across genome for the highest-breadth (pool, virus) hits ----
bg_files <- list.files("results/coverage", pattern = "\\.bedgraph$", full.names = TRUE)
if (length(bg_files) > 0L && nrow(summary) > 0L) {
  read_bg <- function(f) {
    dt <- fread(f, header = FALSE, col.names = c("rname", "start", "end", "depth"))
    dt[, pool_id := sub("\\.bedgraph$", "", basename(f))]
    dt
  }
  bg <- rbindlist(lapply(bg_files, read_bg), use.names = TRUE)
  bg <- merge(bg, refs, by.x = "rname", by.y = "accession", all.x = TRUE)

  # Pick up to 6 best-covered (pool, accession) combinations to visualise.
  top <- head(summary[order(-breadth_pct, -n_reads)], 6L)
  top[, key := paste(pool_id, accession)]
  bg[, key := paste(pool_id, rname)]
  bgp <- bg[key %in% top$key & !is.na(virus)]

  if (nrow(bgp) > 0L) {
    bgp[, facet := paste0(pool_id, "\n", virus,
                          ifelse(is.na(segment) | segment == "NA", "", paste0(" ", segment)))]
    bgp[, mid := (start + end) / 2]
    pc <- ggplot(bgp, aes(mid, depth)) +
      geom_step(direction = "hv", linewidth = 0.3) +
      facet_wrap(~ facet, scales = "free", ncol = 2) +
      labs(title = "Depth across genome (top-covered pool/virus pairs)",
           x = "genome position (bp)", y = "depth") +
      theme_minimal(base_size = 10)
    ggsave("results/figures/coverage_profiles.png", pc, width = 10, height = 7, dpi = 150)
    cat("Wrote results/figures/coverage_profiles.png\n")
  }
}
