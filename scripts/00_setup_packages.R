# =============================================================
# 00. Package setup (run once per computer)
# Installs the Bioconductor packages used in this series.
# =============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("dada2", "ShortRead", "Biostrings"),
                     update = FALSE, ask = FALSE)

# Check they load
library(dada2);      packageVersion("dada2")
library(ShortRead);  packageVersion("ShortRead")