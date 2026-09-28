# =============================================================
# 06. Taxonomic classification (SILVA 138.2, DADA2 naive Bayes)
# Reference: QIIME 2 Moving Pictures tutorial,
#   "qiime feature-classifier classify-sklearn"
# Database: Callahan B. (2024) Silva taxonomic training data formatted
#   for DADA2 (v138.2). Zenodo. doi:10.5281/zenodo.14169026
# =============================================================

suppressPackageStartupMessages(library(dada2))
set.seed(100)
options(timeout = 1800)   # large downloads

# ---- 1. Download SILVA reference files (only once) ----
ref_dir <- "data/reference"
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
zenodo  <- "https://zenodo.org/records/14169026/files/"

refs <- c(
  genus   = "silva_nr99_v138.2_toGenus_trainset.fa.gz",
  species = "silva_v138.2_assignSpecies.fa.gz"
)
md5_expected <- c(
  genus   = "1764e2a36b4500ccb1c7d5261948a414",
  species = "62fd939bd1aa9832d2089847f9f0b87d"
)

for (r in names(refs)) {
  dest <- file.path(ref_dir, refs[[r]])
  if (!file.exists(dest)) {
    download.file(paste0(zenodo, refs[[r]], "?download=1"), dest, mode = "wb")
  }
  ok <- unname(tools::md5sum(dest)) == md5_expected[[r]]
  cat(refs[[r]], "- checksum OK:", ok, "\n")
  if (!ok) stop("Checksum mismatch: delete ", dest, " and run again.")
}

# ---- 2. Load ASVs from Episode 05 ----
seqtab_nochim <- readRDS("data/processed/seqtab_nochim.rds")
asv_seqs <- colnames(seqtab_nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))   # same IDs as Episode 05
cat("ASVs to classify:", length(asv_seqs), "\n")

# ---- 3. Assign taxonomy (Kingdom to Genus) ----
taxa <- assignTaxonomy(seqtab_nochim,
                       file.path(ref_dir, refs[["genus"]]),
                       minBoot = 50, multithread = FALSE, verbose = TRUE)

# ---- 4. Add species (exact matches only) ----
taxa <- addSpecies(taxa, file.path(ref_dir, refs[["species"]]))

# ---- 5. Summaries ----
# Percent of ASVs classified at each rank
pct_classified <- round(100 * colMeans(!is.na(taxa)), 1)
print(pct_classified)

# Top phyla by number of ASVs
print(sort(table(taxa[, "Phylum"], useNA = "ifany"), decreasing = TRUE)[1:10])

# Non-bacterial / organelle sequences (to remove in Episode 07)
cat("Kingdoms:\n"); print(table(taxa[, "Kingdom"], useNA = "ifany"))
cat("Chloroplast ASVs: ", sum(taxa[, "Order"]  == "Chloroplast",  na.rm = TRUE), "\n")
cat("Mitochondria ASVs:", sum(taxa[, "Family"] == "Mitochondria", na.rm = TRUE), "\n")

# ---- 6. Save ----
taxa_table <- data.frame(ASV = asv_ids, taxa, Sequence = asv_seqs,
                         row.names = NULL)
write.csv(taxa_table, "results/tables/06_taxonomy.csv", row.names = FALSE)

rownames(taxa) <- asv_seqs          # keep sequence names for phyloseq later
saveRDS(taxa, "data/processed/taxonomy.rds")