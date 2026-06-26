#!/usr/bin/env Rscript
# Track B: embed reads as canonical tetranucleotide (k=4) frequency vectors,
# reduce with PCA + UMAP, cluster (k-means / DBSCAN / hierarchical), and compare
# the clusters to the minimap2 taxonomy labels.
# Run from the repository root (uses data/raw FASTAs and results/tables/read_hits.tsv).
#
#   Rscript R/40_kmer_clustering.R
#
# Outputs:
#   results/tables/cluster_vs_taxonomy.tsv
#   results/tables/cluster_metrics.tsv
#   results/figures/kmer_umap_taxonomy.png
#   results/figures/kmer_umap_clusters.png

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(uwot)
  library(dbscan)
  library(ggplot2)
})

set.seed(1)
MIN_LEN <- 500L     # drop very short reads: tetranucleotide spectra too noisy
CAP     <- 300L     # cap reads per pool for speed / balance
K_PCS   <- 10L

stopifnot(dir.exists("data"), dir.exists("results"))

fastas <- list.files("data/raw", pattern = "\\.fasta$", full.names = TRUE)
if (length(fastas) == 0L) stop("No FASTAs in data/raw - run scripts/00_fetch_pools.sh first.")

# --- load + subsample reads -----------------------------------------------------
load_pool <- function(f) {
  pool_id <- sub("\\.fasta$", "", basename(f))
  seqs <- readDNAStringSet(f)
  seqs <- seqs[width(seqs) >= MIN_LEN]
  if (length(seqs) == 0L) return(NULL)
  if (length(seqs) > CAP) seqs <- seqs[sample.int(length(seqs), CAP)]
  names(seqs) <- sub("\\s.*$", "", names(seqs))   # keep read id (first token)
  list(seqs = seqs,
       meta = data.table(pool_id = pool_id, read_id = names(seqs)))
}

loaded <- Filter(Negate(is.null), lapply(fastas, load_pool))
if (length(loaded) == 0L) stop("No reads passed the length filter (>=", MIN_LEN, " bp).")

seqs <- do.call(c, lapply(loaded, `[[`, "seqs"))
meta <- rbindlist(lapply(loaded, `[[`, "meta"))

# --- canonical tetranucleotide frequency embedding ------------------------------
of <- oligonucleotideFrequency(seqs, width = 4, as.prob = TRUE)
# Fold each 4-mer with its reverse complement (strand-independent spectrum).
kmers <- colnames(of)
rc    <- as.character(reverseComplement(DNAStringSet(kmers)))
canon <- pmin(kmers, rc)
groups <- split(seq_along(kmers), canon)
emb <- vapply(groups, function(idx) rowSums(of[, idx, drop = FALSE]),
              numeric(nrow(of)))

# Drop reads with no countable k-mers, then zero-variance columns.
ok <- rowSums(emb) > 0 & is.finite(rowSums(emb))
emb  <- emb[ok, , drop = FALSE]
meta <- meta[ok]
emb <- emb[, apply(emb, 2, sd) > 0, drop = FALSE]

# --- taxonomy labels from minimap2 (merge before any position-based assignment) -
# sort = FALSE keeps meta in the same row order as `emb`/`pcs`, so cluster vectors
# (assigned by position below) stay aligned.
if (file.exists("results/tables/read_hits.tsv")) {
  hits <- fread("results/tables/read_hits.tsv")[, .(pool_id, read_id = qname, virus)]
  meta <- merge(meta, hits, by = c("pool_id", "read_id"), all.x = TRUE, sort = FALSE)
}
if (!"virus" %in% names(meta)) meta[, virus := NA_character_]
meta[is.na(virus), virus := "unassigned"]

# --- PCA + UMAP -----------------------------------------------------------------
pca <- prcomp(emb, center = TRUE, scale. = TRUE)
npc <- min(K_PCS, ncol(pca$x))
pcs <- pca$x[, seq_len(npc), drop = FALSE]
n_neighbors <- min(15L, nrow(pcs) - 1L)
set.seed(1)
um  <- umap(pcs, n_neighbors = n_neighbors, min_dist = 0.1, n_threads = 1)
meta[, `:=`(umap1 = um[, 1], umap2 = um[, 2])]

# --- clustering -----------------------------------------------------------------
n_lab <- meta[virus != "unassigned", uniqueN(virus)]
k <- max(2L, min(10L, n_lab))

km  <- kmeans(pcs, centers = k, nstart = 10)
meta[, kmeans := factor(km$cluster)]

hc  <- hclust(dist(pcs), method = "ward.D2")
meta[, hclust := factor(cutree(hc, k = k))]

# DBSCAN with a data-driven eps (median 4-NN distance).
eps <- median(kNNdist(pcs, k = 4))
db  <- dbscan(pcs, eps = eps, minPts = 5)
meta[, dbscan := factor(db$cluster)]   # cluster 0 = noise

# --- concordance: adjusted Rand index (labeled reads only) ----------------------
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  ai <- rowSums(tab); bj <- colSums(tab); n <- sum(tab)
  choose2 <- function(x) x * (x - 1) / 2
  index   <- sum(choose2(tab))
  exp_idx <- sum(choose2(ai)) * sum(choose2(bj)) / choose2(n)
  max_idx <- (sum(choose2(ai)) + sum(choose2(bj))) / 2
  if (max_idx == exp_idx) return(NA_real_)
  (index - exp_idx) / (max_idx - exp_idx)
}

lab <- meta[virus != "unassigned"]
metrics <- data.table(
  method = c("kmeans", "hclust", "dbscan"),
  k_or_eps = c(as.character(k), as.character(k), sprintf("eps=%.3f", eps)),
  n_labeled_reads = nrow(lab),
  adjusted_rand_index = c(
    if (nrow(lab) > 1) adjusted_rand(lab$kmeans, lab$virus) else NA_real_,
    if (nrow(lab) > 1) adjusted_rand(lab$hclust, lab$virus) else NA_real_,
    if (nrow(lab) > 1) adjusted_rand(lab$dbscan, lab$virus) else NA_real_
  )
)
fwrite(metrics, "results/tables/cluster_metrics.tsv", sep = "\t")
cat("Wrote results/tables/cluster_metrics.tsv\n"); print(metrics)

# Cross-tab of k-means cluster vs taxonomy.
xtab <- meta[, .N, by = .(kmeans, virus)]
setorder(xtab, kmeans, -N)
fwrite(xtab, "results/tables/cluster_vs_taxonomy.tsv", sep = "\t")
cat("Wrote results/tables/cluster_vs_taxonomy.tsv\n")

# --- figures --------------------------------------------------------------------
base_theme <- theme_minimal(base_size = 11) + theme(legend.position = "right")

p_tax <- ggplot(meta, aes(umap1, umap2, colour = virus)) +
  geom_point(size = 0.8, alpha = 0.8) +
  labs(title = "k-mer UMAP coloured by minimap2 taxonomy", colour = "virus") +
  base_theme
ggsave("results/figures/kmer_umap_taxonomy.png", p_tax, width = 9, height = 6, dpi = 150)

p_clu <- ggplot(meta, aes(umap1, umap2, colour = kmeans)) +
  geom_point(size = 0.8, alpha = 0.8) +
  labs(title = sprintf("k-mer UMAP coloured by k-means cluster (k=%d)", k),
       colour = "cluster") +
  base_theme
ggsave("results/figures/kmer_umap_clusters.png", p_clu, width = 9, height = 6, dpi = 150)
cat("Wrote UMAP figures to results/figures/\n")
