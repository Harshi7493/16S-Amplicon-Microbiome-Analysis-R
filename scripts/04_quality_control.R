# =============================================================
# 04. Quality control: quality profiles + filtering/trimming
# Reference: QIIME 2 Moving Pictures tutorial,
#   "demux summarize" and the DADA2 settings (trunc-len 120)
# =============================================================

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(ggplot2)
})

# ---- 1. Demultiplexed files from Episode 03 ----
demux_files <- list.files("data/processed/demux", pattern = "\\.fastq\\.gz$",
                          full.names = TRUE)
sample_ids  <- sub("\\.fastq\\.gz$", "", basename(demux_files))
cat("Samples found:", length(demux_files), "\n")

metadata <- read.delim("data/sample-metadata.tsv", comment.char = "#",
                       check.names = FALSE)

# ---- 2. Quality profiles (raw reads) ----
# One example sample per body site
example_ids <- tapply(metadata$`sample-id`, metadata$`body-site`, `[`, 1)
example_files <- demux_files[match(example_ids, sample_ids)]

p_examples <- plotQualityProfile(example_files)
print(p_examples)
ggsave("results/figures/04_quality_examples_raw.png", p_examples,
       width = 10, height = 7, dpi = 300)

# All samples combined
p_all <- plotQualityProfile(demux_files, aggregate = TRUE)
print(p_all)
ggsave("results/figures/04_quality_aggregate_raw.png", p_all,
       width = 8, height = 5, dpi = 300)

# ---- 3. Where are the N (ambiguous) bases? ----
# In Episode 03 we saw an N at position 5 in the first reads.
raw <- readFastq("data/raw/sequences.fastq.gz")
n_by_cycle <- alphabetByCycle(sread(raw))["N", ] / length(raw) * 100
cat("Percent of reads with N, positions 1-10:\n")
print(round(n_by_cycle[1:10], 2))
rm(raw); invisible(gc())

# ---- 4. Filter and trim ----
# Matches the QIIME 2 tutorial: truncate at 120 bp, no left trim.
# maxN = 0   : DADA2 does not allow ambiguous bases
# maxEE = 2  : discard reads with > 2 expected errors
# truncQ = 2 : truncate at the first base with quality <= 2
filt_dir   <- "data/processed/filtered"
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
filt_files <- file.path(filt_dir, basename(demux_files))

filter_out <- filterAndTrim(demux_files, filt_files,
                            truncLen = 120, trimLeft = 0,
                            maxN = 0, maxEE = 2, truncQ = 2,
                            rm.phix = TRUE, compress = TRUE,
                            multithread = FALSE)   # FALSE on Windows

# ---- 5. Filtering summary ----
filter_summary <- data.frame(
  sample_id  = sample_ids,
  reads_in   = filter_out[, "reads.in"],
  reads_out  = filter_out[, "reads.out"],
  pct_kept   = round(100 * filter_out[, "reads.out"] / filter_out[, "reads.in"], 1),
  body_site  = metadata$`body-site`[match(sample_ids, metadata$`sample-id`)]
)
filter_summary <- filter_summary[order(filter_summary$pct_kept), ]

cat("\nTotal reads in: ", sum(filter_summary$reads_in),
    "\nTotal reads out:", sum(filter_summary$reads_out),
    sprintf("\nOverall kept:    %.1f%%\n",
            100 * sum(filter_summary$reads_out) / sum(filter_summary$reads_in)))
print(summary(filter_summary$pct_kept))
print(head(filter_summary, 5))

write.csv(filter_summary, "results/tables/04_filter_summary.csv", row.names = FALSE)

# ---- 6. Quality after filtering ----
p_filt <- plotQualityProfile(filt_files[file.exists(filt_files)], aggregate = TRUE)
print(p_filt)
ggsave("results/figures/04_quality_aggregate_filtered.png", p_filt,
       width = 8, height = 5, dpi = 300)