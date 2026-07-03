# 1. Load Analysis Engines ------------------------------------------------
library(DESeq2)
library(tidyverse)
library(ashr) # Provides empirical Bayes adaptive shrinkage for fold changes

# 2. Data Import & Sample Alignment ---------------------------------------
counts_data <- read.csv("Genes_Quant.csv", row.names = 1, check.names = FALSE)
metadata    <- read.csv("metadata.csv", row.names = 1)

# Enforce matrix columns to match metadata row order exactly
counts_data <- counts_data[, rownames(metadata)]
stopifnot(all(colnames(counts_data) == rownames(metadata)))

# 3. Set Contrast Reference Levels ----------------------------------------
# Defines denominator targets: positive log2FC values imply enrichment in Cooler/2021
metadata$Location <- relevel(factor(metadata$Location), ref = "Hotbox") 
metadata$Year     <- relevel(factor(metadata$Year),     ref = "2019")   

# 4. Construct DESeq2 Dataset Object --------------------------------------
# Additive multi-factor model controlling independently for Year and Location shifts
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_data), # Matrix must consist of raw integer counts
  colData   = metadata,
  design    = ~ Year + Location
)

# 5. Low-Count Sparsity Filtering ----------------------------------------
# Filters for genes with >=1 count in at least 20% of samples AND >=10 total reads
keep <- (rowSums(counts(dds) >= 1) >= ceiling(0.2 * ncol(dds))) & 
  (rowSums(counts(dds)) >= 10)
dds <- dds[keep, ]

# 6. Execute Differential Expression Core ---------------------------------
# Employs 'poscounts' for zero-tolerant size-factor normalization
dds <- DESeq(dds, fitType = "local", sfType = "poscounts")

# Render diagnostic dispersion plot to evaluate fit modeling quality
plotDispEsts(dds)

# 7. Check Modeled Contrast Coefficients ----------------------------------
print(resultsNames(dds))

# 8. Extract & Shrink Year Contrasts (2021 vs 2019) -----------------------
res_year_raw <- results(dds, name = "Year_2021_vs_2019", alpha = 0.05)
res_year_shr <- lfcShrink(dds, coef = "Year_2021_vs_2019", type = "ashr")

summary(res_year_shr)

# 9. Extract & Shrink Location Contrasts (Cooler vs Hotbox) ---------------
res_loc_raw <- results(dds, name = "Location_Cooler_vs_Hotbox", alpha = 0.05)
res_loc_shr <- lfcShrink(dds, coef = "Location_Cooler_vs_Hotbox", type = "ashr")

summary(res_loc_shr)

# 10. File System Export --------------------------------------------------
write.csv(as.data.frame(res_year_shr), file = "DAA_Results_Year_2021_vs_2019.csv")
write.csv(as.data.frame(res_loc_shr),  file = "DAA_Results_Location_Cooler_vs_Hotbox.csv")