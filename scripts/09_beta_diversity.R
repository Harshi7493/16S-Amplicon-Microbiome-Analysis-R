# =============================================================
# 09. Phylogenetic tree, Faith's PD, beta diversity and PCoA
# Reference: QIIME 2 Moving Pictures tutorial,
#   "phylogeny align-to-tree-mafft-fasttree",
#   "diversity core-metrics-phylogenetic" (Jaccard, Bray-Curtis,
#   unweighted/weighted UniFrac, Faith PD, PCoA/Emperor)
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(DECIPHER)
  library(phangorn)
  library(picante)
  library(ggplot2)
})
set.seed(100)

ps_rare <- readRDS("data/processed/ps_rare.rds")   # rarefied to 1103 reads

# ---- 1. Align ASV sequences ----
seqs <- refseq(ps_rare)
cat("ASVs in rarefied table:", length(seqs), "\n")
alignment <- AlignSeqs(seqs, anchor = NA, verbose = FALSE)

# ---- 2. Build tree: neighbour-joining, midpoint-rooted ----
phang <- phyDat(as(alignment, "matrix"), type = "DNA")
dm    <- dist.ml(phang)
tree  <- NJ(dm)
tree$edge.length[tree$edge.length < 0] <- 0     # NJ can give tiny negative lengths
tree  <- midpoint(tree)                          # UniFrac needs a rooted tree

ps_rare <- merge_phyloseq(ps_rare, phy_tree(tree))
ps_rare
ape::write.tree(phy_tree(ps_rare), "results/09_asv_tree.nwk")

# ---- 3. Faith's phylogenetic diversity ----
otu_mat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_mat <- t(otu_mat)
faith <- pd(otu_mat, phy_tree(ps_rare), include.root = TRUE)

md <- data.frame(sample_data(ps_rare))
faith$sample_id <- rownames(faith)
faith$body_site <- md[faith$sample_id, "body_site"]
faith$subject   <- md[faith$sample_id, "subject"]

cat("\nMean Faith's PD by body site:\n")
print(aggregate(PD ~ body_site, data = faith, FUN = function(x) round(mean(x), 2)))
cat("Kruskal-Wallis, Faith PD ~ body site: p =",
    signif(kruskal.test(PD ~ factor(body_site), data = faith)$p.value, 3), "\n")
write.csv(faith, "results/tables/09_faith_pd.csv", row.names = FALSE)

# ---- 4. Distance matrices ----
dists <- list(
  "Jaccard"            = phyloseq::distance(ps_rare, method = "jaccard", binary = TRUE),
  "Bray-Curtis"        = phyloseq::distance(ps_rare, method = "bray"),
  "Unweighted UniFrac" = phyloseq::distance(ps_rare, method = "unifrac"),
  "Weighted UniFrac"   = phyloseq::distance(ps_rare, method = "wunifrac")
)

# ---- 5. PCoA for each metric ----
pcoa_list <- lapply(names(dists), function(n) {
  o  <- ape::pcoa(dists[[n]])
  ve <- round(100 * o$values$Relative_eig[1:2], 1)
  ids <- rownames(o$vectors)
  data.frame(
    sample_id = ids,
    PC1 = o$vectors[, 1], PC2 = o$vectors[, 2],
    metric = n, PC1_pct = ve[1], PC2_pct = ve[2],
    panel  = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)", n, ve[1], ve[2]),
    body_site = md[ids, "body_site"], subject = md[ids, "subject"]
  )
})
pcoa_df <- do.call(rbind, pcoa_list)
pcoa_df$panel <- factor(pcoa_df$panel, levels = unique(pcoa_df$panel))

cat("\nVariance explained by the first two PCoA axes:\n")
print(unique(pcoa_df[, c("metric", "PC1_pct", "PC2_pct")]), row.names = FALSE)

p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, colour = body_site, shape = subject)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ panel, scales = "free") +
  labs(colour = "Body site", shape = "Subject",
       title = "Beta diversity: PCoA (rarefied to 1103 reads)") +
  theme_bw()
print(p_pcoa)
ggsave("results/figures/09_pcoa_4metrics.png", p_pcoa,
       width = 11, height = 8, dpi = 300)

# ---- 6. Save for PERMANOVA (Episode 10) ----
saveRDS(ps_rare, "data/processed/ps_rare_tree.rds")
saveRDS(dists,   "data/processed/beta_distances.rds")