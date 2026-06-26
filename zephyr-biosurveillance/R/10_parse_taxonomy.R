#!/usr/bin/env Rscript
# Parse minimap2 PAF files into one best alignment per read, annotated with the
# reference panel taxonomy. Run from the repository root.
#
#   Rscript R/10_parse_taxonomy.R
#
# Output: results/tables/read_hits.tsv  (one row per mapped read)

suppressPackageStartupMessages({
  library(data.table)
})

stopifnot("Run from repo root (data/ and results/ must exist)" =
            dir.exists("data") && dir.exists("results"))

paf_dir <- "results/paf"
refs    <- fread("data/refs/refs.tsv")
paf_files <- list.files(paf_dir, pattern = "\\.paf$", full.names = TRUE)
if (length(paf_files) == 0L) stop("No PAF files in ", paf_dir, " - run scripts/03_taxonomy_minimap2.sh first.")

# Standard PAF: first 12 columns are fixed; tag columns may follow and are ignored.
paf_cols <- c("qname", "qlen", "qstart", "qend", "strand",
              "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq")

read_paf <- function(f) {
  # PAF has 12 fixed columns plus a variable number of optional tag columns. The
  # varying tag count trips fread's field-count check (it stops early on the first
  # row with an extra tag), so take the 12 fixed columns explicitly via cut to read
  # every alignment line. shQuote handles paths containing spaces.
  dt <- fread(cmd = paste("cut -f1-12", shQuote(f)), header = FALSE, sep = "\t")
  if (nrow(dt) == 0L) return(NULL)
  setnames(dt, paf_cols)
  dt[, pool_id := sub("\\.paf$", "", basename(f))]
  dt[]
}

hits <- rbindlist(lapply(paf_files, read_paf), use.names = TRUE)
if (is.null(hits) || nrow(hits) == 0L) stop("All PAF files were empty (no reads mapped to the panel).")

# Per-alignment identity = residue matches / alignment block length (gap-compressed).
hits[, identity := nmatch / alen]

# Keep the single best alignment per (pool, read): most matched bases, then longest block.
setorder(hits, pool_id, qname, -nmatch, -alen)
best <- hits[, .SD[1L], by = .(pool_id, qname)]

# Annotate with panel taxonomy (tname is the reference accession = first FASTA token).
best <- merge(best, refs, by.x = "tname", by.y = "accession", all.x = TRUE)
best[is.na(virus), virus := "unassigned_in_panel"]

out <- "results/tables/read_hits.tsv"
fwrite(best, out, sep = "\t")
cat(sprintf("Wrote %s (%d mapped reads across %d pools)\n",
            out, nrow(best), uniqueN(best$pool_id)))
