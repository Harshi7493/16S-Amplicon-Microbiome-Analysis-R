# =============================================================
# 03. Importing and inspecting sequence data + demultiplexing
# Reference: QIIME 2 Moving Pictures tutorial, "Demultiplexing"
# (qiime demux emp-single), reproduced here in R.
# =============================================================

suppressPackageStartupMessages({
  library(ShortRead)
  library(Biostrings)
})

# ---- 1. Import raw reads and barcodes ----
reads    <- readFastq("data/raw/sequences.fastq.gz")
barcodes <- readFastq("data/raw/barcodes.fastq.gz")

cat("Sequence reads:", length(reads), "\n")
cat("Barcode reads: ", length(barcodes), "\n")
stopifnot(length(reads) == length(barcodes))   # must pair 1:1

# ---- 2. Inspect ----
print(table(width(reads)))       # read lengths
print(table(width(barcodes)))    # barcode lengths (should be 12)
print(head(sread(reads), 3))     # first few sequences
print(head(id(reads), 3))        # read IDs...
print(head(id(barcodes), 3))     # ...should correspond to these

# ---- 3. Match barcodes to samples ----
metadata <- read.delim("data/sample-metadata.tsv", comment.char = "#",
                       check.names = FALSE)
bc_meta <- setNames(metadata$`barcode-sequence`, metadata$`sample-id`)

bc_reads <- as.character(sread(barcodes))
bc_rc    <- as.character(reverseComplement(DNAStringSet(bc_meta)))
names(bc_rc) <- names(bc_meta)

# Check which orientation matches better
match_fwd <- mean(bc_reads %in% bc_meta)
match_rc  <- mean(bc_reads %in% bc_rc)
cat(sprintf("Barcode match rate - as given: %.1f%%, reverse-complement: %.1f%%\n",
            100 * match_fwd, 100 * match_rc))

lookup <- if (match_rc > match_fwd) bc_rc else bc_meta
sample_of_read <- names(lookup)[match(bc_reads, lookup)]   # NA = no match

# ---- 4. Write one FASTQ file per sample ----
demux_dir <- "data/processed/demux"
dir.create(demux_dir, recursive = TRUE, showWarnings = FALSE)

for (s in names(bc_meta)) {
  out <- file.path(demux_dir, paste0(s, ".fastq.gz"))
  if (file.exists(out)) file.remove(out)
  writeFastq(reads[which(sample_of_read == s)], out, compress = TRUE)
}

# ---- 5. Demultiplexing summary ----
counts <- table(factor(sample_of_read, levels = names(bc_meta)))
reads_per_sample <- data.frame(
  sample_id = names(counts),
  reads     = as.integer(counts),
  body_site = metadata$`body-site`[match(names(counts), metadata$`sample-id`)],
  subject   = metadata$subject[match(names(counts), metadata$`sample-id`)]
)
reads_per_sample <- reads_per_sample[order(reads_per_sample$reads), ]

cat("Unassigned reads:", sum(is.na(sample_of_read)),
    sprintf("(%.1f%%)\n", 100 * mean(is.na(sample_of_read))))
print(summary(reads_per_sample$reads))
print(head(reads_per_sample, 5))    # lowest-depth samples

write.csv(reads_per_sample, "results/tables/03_reads_per_sample.csv",
          row.names = FALSE)