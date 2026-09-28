# =============================================================
# 05. DADA2 denoising and ASV generation
# Reference: QIIME 2 Moving Pictures tutorial,
#   "qiime dada2 denoise-single" (trunc-len 120, consensus chimera removal)
# =============================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(ggplot2)
})
set.seed(100)   # reproducible error learning

# ---- 1. Filtered reads from Episode 04 ----
filt_files <- list.files("data/processed/filtered", pattern = "\\.fastq\\.gz$",
                         full.names = TRUE)
sample_ids <- sub("\\.fastq\\.gz$", "", basename(filt_files))
names(filt_files) <- sample_ids
cat("Filtered samples:", length(filt_files), "\n")

# ---- 2. Learn the error model ----
err <- learnErrors(filt_files, multithread = FALSE)

p_err <- plotErrors(err, nominalQ = TRUE)
print(p_err)
ggsave("results/figures/05_error_model.png", p_err,
       width = 10, height = 8, dpi = 300)

# ---- 3. Denoise each sample (pool = FALSE, as in QIIME 2's default) ----
dada_out <- dada(filt_files, err = err, multithread = FALSE)
dada_out[[1]]   # example summary for one sample

# ---- 4. ASV table ----
seqtab <- makeSequenceTable(dada_out)
cat("ASV table (samples x ASVs):", dim(seqtab), "\n")
print(table(nchar(getSequences(seqtab))))   # should all be 120 bp

# ---- 5. Remove chimeras ----
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = FALSE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab_nochim), "\n")
cat(sprintf("Reads kept after chimera removal: %.1f%%\n",
            100 * sum(seqtab_nochim) / sum(seqtab)))

# ---- 6. Track reads through the pipeline ----
getN <- function(x) sum(getUniques(x))
filter_summary <- read.csv("results/tables/04_filter_summary.csv")

track <- data.frame(
  sample_id    = sample_ids,
  demultiplexed = filter_summary$reads_in[match(sample_ids, filter_summary$sample_id)],
  filtered     = filter_summary$reads_out[match(sample_ids, filter_summary$sample_id)],
  denoised     = sapply(dada_out, getN),
  non_chimeric = rowSums(seqtab_nochim)[sample_ids]
)
track$pct_kept_overall <- round(100 * track$non_chimeric / track$demultiplexed, 1)
track <- track[order(track$non_chimeric), ]
print(head(track, 5))
print(summary(track$non_chimeric))

write.csv(track, "results/tables/05_denoising_stats.csv", row.names = FALSE)

# ---- 7. Save outputs with readable ASV IDs ----
asv_seqs <- getSequences(seqtab_nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))

# Sequences as FASTA
writeXStringSet(DNAStringSet(setNames(asv_seqs, asv_ids)),
                "results/05_asv_sequences.fasta")

# Count table: ASVs as rows, samples as columns (like a QIIME 2 feature table)
asv_table <- t(seqtab_nochim)
rownames(asv_table) <- asv_ids
write.csv(asv_table, "results/tables/05_asv_table.csv")

# R objects for the next episodes
saveRDS(seqtab_nochim, "data/processed/seqtab_nochim.rds")
saveRDS(err,           "data/processed/error_model.rds")