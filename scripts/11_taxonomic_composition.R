# =============================================================
# 11. Taxonomic composition
# Relative abundance bar plots (phylum and genus) per sample and
# per body site, plus summary tables.
# Reference: QIIME 2 Moving Pictures tutorial, "qiime taxa barplot"
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(RColorBrewer)
})

ps_clean <- readRDS("data/processed/ps_clean.rds")   # all 34 samples, unrarefied
suspect  <- c("L2S240", "L3S378")                    # flagged in Episode 09

# ---- 1. Relative abundance, long format ----
ps_rel <- transform_sample_counts(ps_clean, function(x) x / sum(x))
df <- psmelt(ps_rel)

# Readable genus labels for unclassified ASVs
df$Genus_label <- ifelse(!is.na(df$Genus), df$Genus,
                         ifelse(!is.na(df$Family), paste0(df$Family, " (unclassified)"),
                                "Unclassified"))

# Sample labels: order by body site then day; mark suspect samples
df$Sample_label <- ifelse(df$Sample %in% suspect, paste0(df$Sample, "*"), df$Sample)
sample_order <- unique(df[order(df$body_site, df$days_since_experiment_start),
                          "Sample_label"])
df$Sample_label <- factor(df$Sample_label, levels = sample_order)

# ---- 2. Helper: collapse to a rank, keep top N, rest = "Other" ----
collapse_rank <- function(df, rank_col, top_n) {
  per_sample <- aggregate(df$Abundance,
                          by = list(Sample_label = df$Sample_label,
                                    body_site = df$body_site,
                                    Taxon = df[[rank_col]]), FUN = sum)
  names(per_sample)[4] <- "Abundance"
  # Rank taxa by their highest mean abundance in any body site,
  # so taxa dominant in only one site are kept
  site_means <- aggregate(Abundance ~ body_site + Taxon, data = per_sample, FUN = mean)
  max_mean   <- tapply(site_means$Abundance, site_means$Taxon, max)
  top <- names(sort(max_mean, decreasing = TRUE))[1:min(top_n, length(max_mean))]
  per_sample$Taxon <- ifelse(per_sample$Taxon %in% top, per_sample$Taxon, "Other")
  per_sample <- aggregate(Abundance ~ Sample_label + body_site + Taxon,
                          data = per_sample, FUN = sum)
  per_sample$Taxon <- factor(per_sample$Taxon, levels = c(top, "Other"))
  per_sample
}

phylum_df <- collapse_rank(df, "Phylum", 8)
genus_df  <- collapse_rank(df, "Genus_label", 15)

# Colour palette: up to 20 distinct colours, grey for "Other"
make_pal <- function(levels) {
  base <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
  pal  <- setNames(base[seq_len(length(levels) - 1)], head(levels, -1))
  c(pal, Other = "grey75")
}

# ---- 3. Per-sample bar plots ----
bar_plot <- function(d, title) {
  ggplot(d, aes(Sample_label, Abundance, fill = Taxon)) +
    geom_col(width = 0.95) +
    facet_grid(~ body_site, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = make_pal(levels(d$Taxon))) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Relative abundance", fill = NULL, title = title,
         caption = "* suspected cross-contaminated samples (Episode 09)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          legend.text = element_text(size = 8, face = "italic"))
}

p_phylum <- bar_plot(phylum_df, "Phylum-level composition per sample")
print(p_phylum)
ggsave("results/figures/11_barplot_phylum_samples.png", p_phylum,
       width = 12, height = 6, dpi = 300)

p_genus <- bar_plot(genus_df, "Genus-level composition per sample (top 15 genera)")
print(p_genus)
ggsave("results/figures/11_barplot_genus_samples.png", p_genus,
       width = 13, height = 7, dpi = 300)

# ---- 4. Mean composition per body site ----
site_mean <- function(d) {
  m <- aggregate(Abundance ~ body_site + Taxon, data = d, FUN = mean)
  m$Taxon <- factor(m$Taxon, levels = levels(d$Taxon))
  m
}
genus_site <- site_mean(genus_df)

p_site <- ggplot(genus_site, aes(body_site, Abundance, fill = Taxon)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = make_pal(levels(genus_site$Taxon))) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Body site", y = "Mean relative abundance", fill = NULL,
       title = "Mean genus composition by body site") +
  theme_bw() + theme(legend.text = element_text(size = 8, face = "italic"))
print(p_site)
ggsave("results/figures/11_barplot_genus_bodysite.png", p_site,
       width = 8, height = 6, dpi = 300)

# ---- 5. Summary tables (mean % per body site) ----
to_wide <- function(m) {
  w <- round(100 * tapply(m$Abundance, list(m$Taxon, m$body_site), sum), 1)
  w[order(-rowMeans(w)), ]
}
phylum_table <- to_wide(site_mean(phylum_df))
genus_table  <- to_wide(genus_site)

cat("Mean phylum abundance (%) by body site:\n"); print(phylum_table)
cat("\nMean genus abundance (%) by body site (top 15 + Other):\n"); print(genus_table)

write.csv(phylum_table, "results/tables/11_phylum_by_bodysite.csv")
write.csv(genus_table,  "results/tables/11_genus_by_bodysite.csv")