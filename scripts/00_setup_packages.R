# =============================================================
# 00. Package setup (run once per computer)
# Installs all packages used in this series.
# =============================================================

# ---- Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("dada2", "ShortRead", "Biostrings", "phyloseq",
                       "DECIPHER", "ANCOMBC"),
                     update = FALSE, ask = FALSE)

# ---- CRAN packages ----
install.packages(c("phangorn", "picante", "vegan", "ggplot2", "RColorBrewer"))

# ---- Check they load ----
pkgs <- c("dada2", "ShortRead", "Biostrings", "phyloseq", "DECIPHER",
          "ANCOMBC", "phangorn", "picante", "vegan", "ggplot2")
for (p in pkgs) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat(sprintf("%-12s %s\n", p, as.character(packageVersion(p))))
}