#!/usr/bin/env Rscript
# Aggregate per-read hits into per-(pool, virus) calls with high/medium/low confidence,
# combining read support, alignment identity, and genome breadth (from samtools coverage).
# Run from the repository root, after R/10_parse_taxonomy.R and scripts/04_coverage.sh.
#
#   Rscript R/20_confidence_calls.R
#
# Outputs:
#   results/tables/taxonomy_calls.tsv
#   results/figures/taxonomy_heatmap.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

stopifnot(dir.exists("data"), dir.exists("results"))

refs  <- fread("data/refs/refs.tsv")
hits  <- fread("results/tables/read_hits.tsv")

# --- breadth/depth per (pool, virus) from samtools coverage files ---------------
cov_files <- list.files("results/coverage", pattern = "\\.coverage\\.txt$", full.names = TRUE)
cov_by_virus <- NULL
if (length(cov_files) > 0L) {
  read_cov <- function(f) {
    dt <- fread(f)
    setnames(dt, 1, "rname")                 # samtools writes "#rname" as first column
    dt[, pool_id := sub("\\.coverage\\.txt$", "", basename(f))]
    dt
  }
  cov <- rbindlist(lapply(cov_files, read_cov), use.names = TRUE, fill = TRUE)
  cov <- merge(cov, refs, by.x = "rname", by.y = "accession", all.x = TRUE)
  # Per virus: best breadth across its references, and the matching mean depth.
  cov_by_virus <- cov[!is.na(virus),
    .(breadth_pct = max(coverage, na.rm = TRUE),
      mean_depth  = max(meandepth, na.rm = TRUE)),
    by = .(pool_id, virus)]
}

# --- read support + identity per (pool, virus) ----------------------------------
mapped <- hits[virus != "unassigned_in_panel"]
calls <- mapped[, .(
  n_reads        = .N,
  median_identity = round(median(identity), 3),
  total_aln_bp   = sum(alen),
  n_refs_hit     = uniqueN(tname)
), by = .(pool_id, virus)]

# Family, for ambiguity flagging among close relatives.
fam <- unique(refs[, .(virus, family)])
calls <- merge(calls, fam, by = "virus", all.x = TRUE)

if (!is.null(cov_by_virus)) {
  calls <- merge(calls, cov_by_virus, by = c("pool_id", "virus"), all.x = TRUE)
} else {
  calls[, `:=`(breadth_pct = NA_real_, mean_depth = NA_real_)]
}
calls[is.na(breadth_pct), breadth_pct := 0]
calls[is.na(mean_depth),  mean_depth  := 0]

# --- confidence tiers (heuristic; see memo/README) ------------------------------
classify <- function(n_reads, median_identity, breadth_pct) {
  if (n_reads >= 50 && median_identity >= 0.90 && breadth_pct >= 25) return("high")
  if (n_reads >= 5  && median_identity >= 0.80 && breadth_pct >= 5)  return("medium")
  "low"
}
calls[, confidence := mapply(classify, n_reads, median_identity, breadth_pct)]
calls[, confidence := factor(confidence, levels = c("high", "medium", "low"))]

# Flag pools where multiple members of the same viral family are called (possible
# cross-mapping / ambiguity among close relatives, e.g. seasonal coronaviruses).
calls[, n_family_in_pool := uniqueN(virus), by = .(pool_id, family)]
calls[, ambiguity_flag := n_family_in_pool > 1]

setorder(calls, pool_id, confidence, -n_reads)
out <- "results/tables/taxonomy_calls.tsv"
fwrite(calls, out, sep = "\t")
cat(sprintf("Wrote %s (%d pool-virus calls)\n", out, nrow(calls)))

# --- heatmap: virus x pool, fill = log10 reads, outlined by confidence ----------
calls[, log10_reads := log10(n_reads)]
p <- ggplot(calls, aes(x = pool_id, y = virus, fill = log10_reads)) +
  geom_tile(aes(colour = confidence), linewidth = 0.6, width = 0.95, height = 0.95) +
  scale_fill_viridis_c(option = "C", name = "log10(reads)") +
  scale_colour_manual(values = c(high = "black", medium = "grey40", low = "grey80"),
                      name = "confidence", drop = FALSE) +
  labs(title = "Viral reads per pool (panel references)",
       x = "pool", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

fig <- "results/figures/taxonomy_heatmap.png"
ggsave(fig, p, width = 10, height = 6, dpi = 150)
cat(sprintf("Wrote %s\n", fig))
