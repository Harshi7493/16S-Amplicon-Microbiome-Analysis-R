# =============================================================
# 07. Feature table exploration
# Builds a phyloseq object (ASV table + taxonomy + metadata),
# removes non-target sequences, and summarises the table.
# Reference: QIIME 2 Moving Pictures tutorial,
#   "feature-table summarize" / "feature-table tabulate-seqs"
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(ggplot2)
})

# ---- 1. Load inputs ----
seqtab_nochim <- readRDS("data/processed/seqtab_nochim.rds")
taxa          <- readRDS("data/processed/taxonomy.rds")

metadata <- read.delim("data/sample-metadata.tsv", comment.char = "#",
                       check.names = FALSE)
names(metadata) <- gsub("-", "_", names(metadata))   # e.g. body-site -> body_site
rownames(metadata) <- metadata$sample_id

# ---- 2. Build the phyloseq object ----
ps <- phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

# Store sequences separately and use short ASV names
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq_len(ntaxa(ps)))
ps

# ---- 3. Remove non-target sequences ----
tax <- as.data.frame(tax_table(ps))
is_mito   <- tax$Family  %in% "Mitochondria"
is_chloro <- tax$Order   %in% "Chloroplast"
is_euk    <- tax$Kingdom %in% "Eukaryota"
no_phylum <- is.na(tax$Phylum)
to_remove <- is_mito | is_chloro | is_euk | no_phylum

removal_summary <- data.frame(
  category = c("Mitochondria", "Chloroplast", "Eukaryota", "No phylum", "Total removed"),
  ASVs  = c(sum(is_mito), sum(is_chloro), sum(is_euk), sum(no_phylum), sum(to_remove)),
  reads = c(sum(taxa_sums(ps)[is_mito]), sum(taxa_sums(ps)[is_chloro]),
            sum(taxa_sums(ps)[is_euk]),  sum(taxa_sums(ps)[no_phylum]),
            sum(taxa_sums(ps)[to_remove]))
)
print(removal_summary)

ps_clean <- prune_taxa(!to_remove, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
ps_clean

# Share of each sample's reads that was removed
pct_removed <- round(100 * (1 - sample_sums(ps_clean) / sample_sums(ps)), 1)
cat("\nPercent of reads removed per sample:\n")
print(summary(pct_removed))
print(tapply(pct_removed, sample_data(ps)$body_site, mean))   # by body site

# ---- 4. Per-sample summary ----
sample_summary <- data.frame(
  sample_id    = sample_names(ps_clean),
  body_site    = sample_data(ps_clean)$body_site,
  subject      = sample_data(ps_clean)$subject,
  reads        = sample_sums(ps_clean),
  observed_asvs = colSums(t(otu_table(ps_clean)) > 0)
)
sample_summary <- sample_summary[order(sample_summary$reads), ]
cat("\nReads per sample:\n");  print(summary(sample_summary$reads))
cat("\nMean reads and ASVs by body site:\n")
print(aggregate(cbind(reads, observed_asvs) ~ body_site, data = sample_summary, FUN = mean))

write.csv(sample_summary, "results/tables/07_sample_summary.csv", row.names = FALSE)

# ---- 5. Per-ASV summary: abundance and prevalence ----
prevalence <- colSums(otu_table(ps_clean) > 0)   # number of samples containing each ASV
asv_summary <- data.frame(
  ASV        = taxa_names(ps_clean),
  total_reads = taxa_sums(ps_clean),
  prevalence = prevalence,
  as.data.frame(tax_table(ps_clean))[, c("Phylum", "Family", "Genus")]
)
asv_summary <- asv_summary[order(-asv_summary$total_reads), ]

cat("\nASVs found in only one sample:", sum(prevalence == 1),
    sprintf("(%.1f%% of ASVs)\n", 100 * mean(prevalence == 1)))
cat("\nTop 10 ASVs by total reads:\n")
print(head(asv_summary, 10), row.names = FALSE)

write.csv(asv_summary, "results/tables/07_asv_summary.csv", row.names = FALSE)

# ---- 6. Plot: sequencing depth per sample ----
p_depth <- ggplot(sample_summary, aes(x = body_site, y = reads, colour = subject)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
  geom_hline(yintercept = 1103, linetype = "dashed") +
  annotate("text", x = 0.6, y = 1103, label = "1,103 (tutorial depth)",
           vjust = -0.5, hjust = 0, size = 3.5) +
  labs(x = "Body site", y = "Reads per sample (after filtering)",
       colour = "Subject",
       title = "Sequencing depth per sample") +
  theme_bw()
print(p_depth)
ggsave("results/figures/07_depth_per_sample.png", p_depth,
       width = 7, height = 5, dpi = 300)

# ---- 7. Save the cleaned object for the diversity episodes ----
saveRDS(ps_clean, "data/processed/ps_clean.rds")