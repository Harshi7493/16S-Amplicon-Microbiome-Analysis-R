# =============================================================
# 08. Alpha diversity
# Rarefaction curves, rarefying to an even depth, richness/diversity
# metrics, and statistical tests across body sites.
# Reference: QIIME 2 Moving Pictures tutorial,
#   "core-metrics-phylogenetic --p-sampling-depth 1103",
#   "diversity alpha-group-significance", "alpha-correlation"
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})

ps_clean <- readRDS("data/processed/ps_clean.rds")
depth_main <- 1103   # QIIME 2 tutorial sampling depth
depth_sens <- 824    # lowest sample: sensitivity check

# ---- 1. Rarefaction curves ----
otu <- as(otu_table(ps_clean), "matrix")   # samples x ASVs
rc  <- rarecurve(otu, step = 50, tidy = TRUE)
names(rc) <- c("sample_id", "reads", "ASVs")
rc$body_site <- sample_data(ps_clean)$body_site[match(rc$sample_id, sample_names(ps_clean))]

p_rare <- ggplot(rc, aes(reads, ASVs, group = sample_id, colour = body_site)) +
  geom_line(alpha = 0.8) +
  geom_vline(xintercept = depth_main, linetype = "dashed") +
  labs(x = "Reads sampled", y = "Observed ASVs", colour = "Body site",
       title = "Rarefaction curves") +
  theme_bw()
print(p_rare)
ggsave("results/figures/08_rarefaction_curves.png", p_rare,
       width = 8, height = 5, dpi = 300)

# ---- 2. Rarefy to even depth ----
ps_rare <- rarefy_even_depth(ps_clean, sample.size = depth_main,
                             rngseed = 100, replace = FALSE, trimOTUs = TRUE)
cat("Samples kept at", depth_main, "reads:", nsamples(ps_rare), "\n")
cat("Dropped:", setdiff(sample_names(ps_clean), sample_names(ps_rare)), "\n")

# ---- 3. Alpha diversity metrics ----
alpha_metrics <- function(ps) {
  a <- estimate_richness(ps, measures = c("Observed", "Shannon"))
  a$Pielou    <- a$Shannon / log(a$Observed)        # evenness
  a$sample_id <- rownames(a)
  md <- data.frame(sample_data(ps))
  cbind(a, md[a$sample_id, c("body_site", "subject",
                             "reported_antibiotic_usage",
                             "days_since_experiment_start")])
}
alpha <- alpha_metrics(ps_rare)
alpha$body_site <- factor(alpha$body_site)

cat("\nMean alpha diversity by body site:\n")
print(aggregate(cbind(Observed, Shannon, Pielou) ~ body_site, data = alpha,
                FUN = function(x) round(mean(x), 2)))

write.csv(alpha, "results/tables/08_alpha_diversity.csv", row.names = FALSE)

# ---- 4. Statistical tests ----
# Kruskal-Wallis across body sites (as in alpha-group-significance)
kw <- data.frame(
  metric = c("Observed", "Shannon", "Pielou"),
  p_body_site = sapply(c("Observed", "Shannon", "Pielou"), function(m)
    kruskal.test(alpha[[m]] ~ alpha$body_site)$p.value),
  p_subject = sapply(c("Observed", "Shannon", "Pielou"), function(m)
    kruskal.test(alpha[[m]] ~ factor(alpha$subject))$p.value)
)
cat("\nKruskal-Wallis p-values:\n"); print(kw, digits = 3)

# Pairwise body-site comparisons for Shannon (BH-adjusted)
cat("\nPairwise Wilcoxon (Shannon), BH-adjusted p-values:\n")
print(pairwise.wilcox.test(alpha$Shannon, alpha$body_site,
                           p.adjust.method = "BH", exact = FALSE))

# Correlation with time (as in alpha-correlation)
cat("\nSpearman correlation, Shannon vs days since start:\n")
print(cor.test(alpha$Shannon, alpha$days_since_experiment_start,
               method = "spearman", exact = FALSE))

write.csv(kw, "results/tables/08_alpha_kruskal.csv", row.names = FALSE)

# ---- 5. Plot ----
alpha_long <- rbind(
  data.frame(alpha[, c("sample_id", "body_site", "subject")], metric = "Observed ASVs", value = alpha$Observed),
  data.frame(alpha[, c("sample_id", "body_site", "subject")], metric = "Shannon",       value = alpha$Shannon),
  data.frame(alpha[, c("sample_id", "body_site", "subject")], metric = "Pielou evenness", value = alpha$Pielou)
)
alpha_long$metric <- factor(alpha_long$metric,
                            levels = c("Observed ASVs", "Shannon", "Pielou evenness"))

p_alpha <- ggplot(alpha_long, aes(body_site, value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(aes(colour = subject), width = 0.15, size = 2.5, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Body site", y = NULL, colour = "Subject",
       title = paste0("Alpha diversity (rarefied to ", depth_main, " reads)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_alpha)
ggsave("results/figures/08_alpha_diversity.png", p_alpha,
       width = 10, height = 4.5, dpi = 300)

# ---- 6. Sensitivity check at the lowest depth (all samples) ----
ps_rare_sens <- rarefy_even_depth(ps_clean, sample.size = depth_sens,
                                  rngseed = 100, replace = FALSE, trimOTUs = TRUE)
alpha_sens <- alpha_metrics(ps_rare_sens)
cat("\nSensitivity check at", depth_sens, "reads (", nsamples(ps_rare_sens), "samples ):\n")
print(sapply(c("Observed", "Shannon", "Pielou"), function(m)
  kruskal.test(alpha_sens[[m]] ~ factor(alpha_sens$body_site))$p.value), digits = 3)

# ---- 7. Save rarefied object for beta diversity ----
saveRDS(ps_rare, "data/processed/ps_rare.rds")