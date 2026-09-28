# =============================================================
# 13. Publication-quality multi-panel figure
# Combines alpha diversity, PCoA, composition and differential
# abundance into Figure 1 (PNG for README, PDF for print).
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ape)
  library(RColorBrewer)
})

flagged <- c("L2S240", "L3S378", "L3S242")   # day-0 palm samples (Episodes 09-11)

# Colour-blind-friendly (Okabe-Ito) colours, consistent across panels
site_cols <- c("gut" = "#E69F00", "left palm" = "#56B4E9",
               "right palm" = "#0072B2", "tongue" = "#CC79A7")

theme_pub <- theme_classic(base_size = 9) +
  theme(plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8),
        legend.title  = element_text(size = 8),
        legend.text   = element_text(size = 7),
        legend.key.size = unit(3.5, "mm"))

# ---- Panel A: Shannon diversity ----
alpha <- read.csv("results/tables/08_alpha_diversity.csv")
alpha$flagged <- alpha$sample_id %in% flagged
kw <- read.csv("results/tables/08_alpha_kruskal.csv")
p_shannon <- kw$p_body_site[kw$metric == "Shannon"]

pA <- ggplot(alpha, aes(body_site, Shannon)) +
  geom_boxplot(aes(fill = body_site), alpha = 0.45, width = 0.6,
               outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(aes(colour = body_site, shape = flagged), width = 0.12, size = 1.8) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_colour_manual(values = site_cols, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  labs(x = NULL, y = "Shannon diversity", title = "Alpha diversity",
       subtitle = sprintf("Kruskal-Wallis p = %.3f (rarefied to 1,103 reads)", p_shannon)) +
  theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ---- Panel B: PCoA, weighted UniFrac ----
dists <- readRDS("data/processed/beta_distances.rds")
pc_res <- pcoa(dists[["Weighted UniFrac"]])
ve <- 100 * pc_res$values$Relative_eig[1:2]
pc <- data.frame(sample_id = rownames(pc_res$vectors),
                 PC1 = pc_res$vectors[, 1], PC2 = pc_res$vectors[, 2])
pc <- merge(pc, alpha[, c("sample_id", "body_site", "subject")], by = "sample_id")
pc$flagged <- pc$sample_id %in% flagged

perm <- read.csv("results/tables/10_permanova_main.csv")
wu <- perm[perm$metric == "Weighted UniFrac" & perm$term == "body_site", ]

pB <- ggplot(pc, aes(PC1, PC2)) +
  geom_point(aes(colour = body_site, shape = subject), size = 2.2, alpha = 0.9) +
  geom_point(data = pc[pc$flagged, ], shape = 21, size = 4.2,
             colour = "black", fill = NA, stroke = 0.5) +
  scale_colour_manual(values = site_cols) +
  scale_shape_manual(values = c("subject-1" = 16, "subject-2" = 17)) +
  labs(x = sprintf("PCoA1 (%.1f%%)", ve[1]), y = sprintf("PCoA2 (%.1f%%)", ve[2]),
       colour = "Body site", shape = "Subject",
       title = "Beta diversity (weighted UniFrac)",
       subtitle = sprintf("PERMANOVA body site R\u00b2 = %.2f, p = %.3f; circled = flagged samples",
                          wu$R2, wu$p)) +
  theme_pub

# ---- Panel C: mean genus composition ----
g <- as.matrix(read.csv("results/tables/11_genus_by_bodysite.csv",
                        row.names = 1, check.names = FALSE))
g_other <- g["Other", ]
g_taxa  <- g[rownames(g) != "Other", ]
top10   <- names(sort(apply(g_taxa, 1, max), decreasing = TRUE))[1:10]
gm <- rbind(g_taxa[top10, ],
            Other = colSums(g_taxa[setdiff(rownames(g_taxa), top10), , drop = FALSE]) + g_other)

comp <- data.frame(Genus = rep(rownames(gm), times = ncol(gm)),
                   body_site = rep(colnames(gm), each = nrow(gm)),
                   pct = as.vector(gm))
comp$Genus <- factor(comp$Genus, levels = rownames(gm))
genus_cols <- setNames(c(brewer.pal(10, "Paired"), "grey80"), rownames(gm))

pC <- ggplot(comp, aes(body_site, pct, fill = Genus)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = genus_cols) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Genus",
       title = "Taxonomic composition", subtitle = "Top 10 genera; all 34 samples") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.text = element_text(size = 7, face = "italic"))

# ---- Panel D: differential abundance, skin vs tongue ----
da <- read.csv("results/tables/12_ancombc_skin_vs_tongue.csv")
da$taxon <- sub("^Genus:", "", da$taxon)
s <- da[da$significant %in% TRUE, ]
s$higher_in <- ifelse(s$lfc > 0, "skin", "tongue")
s$taxon <- factor(s$taxon, levels = s$taxon[order(s$lfc)])

pD <- ggplot(s, aes(lfc, taxon, fill = higher_in)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = lfc - se, xmax = lfc + se), width = 0.25, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(skin = "#0072B2", tongue = "#CC79A7")) +
  labs(x = "Log fold change (\u00b1 SE)", y = NULL, fill = "Higher in",
       title = "Differential abundance: skin vs tongue",
       subtitle = "ANCOM-BC2, Holm q < 0.05; flagged samples excluded") +
  theme_pub + theme(axis.text.y = element_text(face = "italic"))

# ---- Assemble and export ----
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

print(fig)
ggsave("results/figures/13_figure1.png", fig,
       width = 180, height = 170, units = "mm", dpi = 300, bg = "white")
ggsave("results/figures/13_figure1.pdf", fig,
       width = 180, height = 170, units = "mm")
cat("Saved: results/figures/13_figure1.png and .pdf\n")