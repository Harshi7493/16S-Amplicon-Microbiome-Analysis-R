# =============================================================
# 02. Moving Pictures dataset & metadata
# Downloads the QIIME 2 (2024.10) Moving Pictures tutorial data
# and explores the sample metadata.
# =============================================================

options(timeout = 600)  # allow time for larger downloads

# ---- Folders (data/raw is git-ignored, so create it here) ----
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

base_url <- "https://data.qiime2.org/2024.10/tutorials/moving-pictures/"

# ---- Download metadata (small: tracked in Git) ----
download.file(paste0(base_url, "sample_metadata.tsv"),
              destfile = "data/sample-metadata.tsv")

# ---- Download raw sequences (large: git-ignored) ----
download.file(paste0(base_url, "emp-single-end-sequences/barcodes.fastq.gz"),
              destfile = "data/raw/barcodes.fastq.gz", mode = "wb")
download.file(paste0(base_url, "emp-single-end-sequences/sequences.fastq.gz"),
              destfile = "data/raw/sequences.fastq.gz", mode = "wb")

file.info(list.files("data/raw", full.names = TRUE))[, "size", drop = FALSE]

# ---- Read metadata ----
# Row 2 of the file ("#q2:types") holds QIIME 2 column types,
# so comment.char = "#" skips it.
metadata <- read.delim("data/sample-metadata.tsv", comment.char = "#",
                       check.names = FALSE)

# ---- Explore the study design ----
dim(metadata)                                   # samples x variables
str(metadata)
table(metadata$`body-site`)                     # 4 body sites
table(metadata$subject)                         # 2 subjects
table(metadata$`reported-antibiotic-usage`)     # antibiotic exposure
table(metadata$`body-site`, metadata$subject)   # design balance