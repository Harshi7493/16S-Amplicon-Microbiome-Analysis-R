# =============================================================
# 10. PERMANOVA, dispersion tests and sensitivity analysis
# Reference: QIIME 2 Moving Pictures tutorial,
#   "diversity beta-group-significance" (pairwise PERMANOVA)
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})
set.seed(100)

ps_rare <- readRDS("data/processed/ps_rare_tree.rds")
dists   <- readRDS("data/processed/beta_distances.rds")
md      <- data.frame(sample_data(ps_rare))

suspect <- c("L2S240", "L3S378")   # day-0 palm samples flagged in Episode 09

# Helper: drop samples from a dist object
drop_samples <- function(d, drop) {
  m <- as.matrix(d); keep <- !rownames(m) %in% drop
  as.dist(m[keep, keep])
}

# ---- 1. PERMANOVA: body site + subject (marginal effects) ----
run_permanova <- function(d) {
  m <- md[labels(d), ]
  r <- adonis2(d ~ body_site + subject, data = m, permutations = 999, by = "margin")
  data.frame(term = rownames(r)[1:2], R2 = round(r$R2[1:2], 3),
             F = round(r$F[1:2], 2), p = r$`Pr(>F)`[1:2])
}

permanova_main <- do.call(rbind, lapply(names(dists), function(n)
  cbind(metric = n, run_permanova(dists[[n]]))))
cat("PERMANOVA (all 31 samples):\n"); print(permanova_main, row.names = FALSE)
write.csv(permanova_main, "results/tables/10_permanova_main.csv", row.names = FALSE)

# ---- 2. Dispersion (betadisper) by body site ----
disp_list <- lapply(names(dists), function(n) {
  d  <- dists[[n]]
  bd <- betadisper(d, md[labels(d), "body_site"])
  pt <- permutest(bd, permutations = 999)
  list(test = data.frame(metric = n, F = round(pt$tab$F[1], 2), p = pt$tab$`Pr(>F)`[1]),
       dist = data.frame(metric = n, body_site = bd$group, dist_to_centroid = bd$distances))
})
disp_tests <- do.call(rbind, lapply(disp_list, `[[`, "test"))
cat("\nDispersion test (betadisper) by body site:\n"); print(disp_tests, row.names = FALSE)
write.csv(disp_tests, "results/tables/10_betadisper.csv", row.names = FALSE)

disp_df <- do.call(rbind, lapply(disp_list, `[[`, "dist"))
disp_df$metric <- factor(disp_df$metric, levels = names(dists))
p_disp <- ggplot(disp_df, aes(body_site, dist_to_centroid)) +
  geom_boxplot(fill = "grey90", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Body site", y = "Distance to group centroid",
       title = "Within-group dispersion (betadisper)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_disp)
ggsave("results/figures/10_dispersion.png", p_disp, width = 10, height = 7, dpi = 300)

# ---- 3. Pairwise PERMANOVA between body sites (BH-adjusted) ----
pairwise_permanova <- function(d) {
  m   <- as.matrix(d)
  grp <- md[rownames(m), "body_site"]
  prs <- combn(sort(unique(grp)), 2)
  out <- do.call(rbind, lapply(seq_len(ncol(prs)), function(i) {
    p    <- prs[, i]; keep <- grp %in% p
    dd   <- as.dist(m[keep, keep])
    r    <- adonis2(dd ~ g, data = data.frame(g = factor(grp[keep])), permutations = 999)
    data.frame(group1 = p[1], group2 = p[2], R2 = round(r$R2[1], 3),
               F = round(r$F[1], 2), p = r$`Pr(>F)`[1])
  }))
  out$p_adj <- round(p.adjust(out$p, method = "BH"), 4)
  out
}

pairwise_uu <- pairwise_permanova(dists[["Unweighted UniFrac"]])
cat("\nPairwise PERMANOVA - Unweighted UniFrac (as in the tutorial):\n")
print(pairwise_uu, row.names = FALSE)
pairwise_bc <- pairwise_permanova(dists[["Bray-Curtis"]])
cat("\nPairwise PERMANOVA - Bray-Curtis:\n")
print(pairwise_bc, row.names = FALSE)

write.csv(cbind(metric = "Unweighted UniFrac", pairwise_uu),
          "results/tables/10_pairwise_unweighted_unifrac.csv", row.names = FALSE)
write.csv(cbind(metric = "Bray-Curtis", pairwise_bc),
          "results/tables/10_pairwise_bray_curtis.csv", row.names = FALSE)

# ---- 4. Stricter test: permute body site only within each subject ----
restricted <- do.call(rbind, lapply(names(dists), function(n) {
  d <- dists[[n]]; m <- md[labels(d), ]
  r <- adonis2(d ~ body_site, data = m,
               permutations = how(nperm = 999, blocks = factor(m$subject)))
  data.frame(metric = n, R2 = round(r$R2[1], 3), F = round(r$F[1], 2), p = r$`Pr(>F)`[1])
}))
cat("\nBody site, permutations restricted within subject:\n")
print(restricted, row.names = FALSE)

# ---- 5. Sensitivity: exclude the two flagged day-0 palm samples ----
permanova_sens <- do.call(rbind, lapply(names(dists), function(n)
  cbind(metric = n, run_permanova(drop_samples(dists[[n]], suspect)))))
cat("\nPERMANOVA without", paste(suspect, collapse = " and "), "(29 samples):\n")
print(permanova_sens, row.names = FALSE)
write.csv(permanova_sens, "results/tables/10_permanova_sensitivity.csv", row.names = FALSE)