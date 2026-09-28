# =============================================================
# 12. Differential abundance with ANCOM-BC2
# Part A: gut samples, subject-1 vs subject-2 (as in the tutorial)
# Part B: skin (palms) vs tongue, adjusted for subject (extension)
# Reference: QIIME 2 Moving Pictures tutorial,
#   "qiime composition ancombc" (gut samples, formula 'subject')
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
})
set.seed(100)

ps_clean <- readRDS("data/processed/ps_clean.rds")

# ---- 1. Exclude the day-0 palm samples flagged as cross-contaminated ----
flagged <- c("L2S240", "L3S378", "L3S242")
ps_da <- prune_samples(!sample_names(ps_clean) %in% flagged, ps_clean)
cat("Samples used:", nsamples(ps_da), "(excluded:", paste(flagged, collapse = ", "), ")\n")

# ---- Helper: tidy one coefficient from ancombc2 output ----
tidy_ancombc <- function(res, var) {
  lfc_col <- grep(paste0("^lfc_", var), names(res), value = TRUE)[1]
  suffix  <- sub("^lfc_", "", lfc_col)
  out <- data.frame(
    taxon     = res$taxon,
    lfc       = res[[paste0("lfc_", suffix)]],
    se        = res[[paste0("se_", suffix)]],
    q         = res[[paste0("q_", suffix)]],
    diff      = res[[paste0("diff_", suffix)]],
    passed_ss = res[[paste0("passed_ss_", suffix)]]
  )
  out$significant <- out$diff & out$passed_ss   # significant AND robust to pseudo-counts
  out[order(out$q), ]
}

# ---- Helper: plot significant log-fold changes ----
plot_lfc <- function(d, title, pos_label, neg_label) {
  s <- d[d$significant, ]
  if (nrow(s) == 0) { message("No significant taxa for: ", title); return(NULL) }
  s$direction <- ifelse(s$lfc > 0, pos_label, neg_label)
  s$taxon <- factor(s$taxon, levels = s$taxon[order(s$lfc)])
  ggplot(s, aes(lfc, taxon, fill = direction)) +
    geom_col() +
    geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), height = 0.3) +
    geom_vline(xintercept = 0) +
    labs(x = "Log fold change (natural log, ± SE)", y = NULL, fill = "Higher in",
         title = title) +
    theme_bw() + theme(axis.text.y = element_text(face = "italic"))
}

# ================= PART A: gut, subject-1 vs subject-2 =================
ps_gut <- subset_samples(ps_da, body_site == "gut")
ps_gut <- prune_taxa(taxa_sums(ps_gut) > 0, ps_gut)
sample_data(ps_gut)$subject <- factor(sample_data(ps_gut)$subject,
                                      levels = c("subject-1", "subject-2"))
cat("\nPart A - gut samples:", nsamples(ps_gut), "\n")

out_gut <- ancombc2(data = ps_gut, tax_level = "Genus",
                    fix_formula = "subject", p_adj_method = "holm",
                    prv_cut = 0.10, lib_cut = 0, group = NULL,
                    struc_zero = FALSE, alpha = 0.05, n_cl = 1, verbose = FALSE)

res_gut <- tidy_ancombc(out_gut$res, "subject")
cat("Genera tested:", nrow(res_gut), "| significant:", sum(res_gut$significant, na.rm = TRUE), "\n")
print(head(res_gut[, c("taxon", "lfc", "se", "q", "significant")], 10), row.names = FALSE)
write.csv(res_gut, "results/tables/12_ancombc_gut_subject.csv", row.names = FALSE)

p_gut <- plot_lfc(res_gut, "Gut genera: subject-2 vs subject-1 (ANCOM-BC2)",
                  "subject-2", "subject-1")
if (!is.null(p_gut)) {
  print(p_gut)
  ggsave("results/figures/12_ancombc_gut_subject.png", p_gut,
         width = 8, height = max(3, 0.35 * sum(res_gut$significant) + 1.5), dpi = 300)
}

# ================= PART B: skin (palms) vs tongue =================
ps_st <- subset_samples(ps_da, body_site != "gut")
ps_st <- prune_taxa(taxa_sums(ps_st) > 0, ps_st)
sample_data(ps_st)$habitat <- factor(
  ifelse(sample_data(ps_st)$body_site == "tongue", "tongue", "skin"),
  levels = c("tongue", "skin"))                 # reference = tongue
sample_data(ps_st)$subject <- factor(sample_data(ps_st)$subject)
cat("\nPart B - samples:", nsamples(ps_st), "\n")
print(table(sample_data(ps_st)$habitat))

out_st <- ancombc2(data = ps_st, tax_level = "Genus",
                   fix_formula = "habitat + subject", p_adj_method = "holm",
                   prv_cut = 0.10, lib_cut = 0, group = NULL,
                   struc_zero = FALSE, alpha = 0.05, n_cl = 1, verbose = FALSE)

res_st <- tidy_ancombc(out_st$res, "habitat")
cat("Genera tested:", nrow(res_st), "| significant:", sum(res_st$significant, na.rm = TRUE), "\n")
cat("\nHigher on skin (top 10):\n")
print(head(res_st[res_st$significant & res_st$lfc > 0,
                  c("taxon", "lfc", "se", "q")], 10), row.names = FALSE)
cat("\nHigher on tongue (top 10):\n")
print(head(res_st[res_st$significant & res_st$lfc < 0,
                  c("taxon", "lfc", "se", "q")], 10), row.names = FALSE)
write.csv(res_st, "results/tables/12_ancombc_skin_vs_tongue.csv", row.names = FALSE)

p_st <- plot_lfc(res_st, "Genera differing between skin and tongue (ANCOM-BC2)",
                 "skin", "tongue")
if (!is.null(p_st)) {
  print(p_st)
  ggsave("results/figures/12_ancombc_skin_vs_tongue.png", p_st,
         width = 8, height = max(3, 0.3 * sum(res_st$significant) + 1.5), dpi = 300)
}