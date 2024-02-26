## Concatenate pseudobulk gene score matrices of different chromosomes

# run with conda environment DESeq2

library(dplyr)
library(data.table)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# normalized or not
# norm <- "_Normalized"
norm <- ""

# 1. Read pseudobulk gene score matrices for the different chromosomes
p <- paste0("Pseudobulk_GeneScoreMatrix", norm, "_chr.*.csv")
files <- list.files(path = paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/"),
            pattern = p,
            full.names = FALSE)
chrom_list <- list()
for (f in files){
    chrom_mat <- fread(paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/",f))
    chrom_list[[f]] <- as.data.frame(t(chrom_mat))
}
chrom_df <- do.call(rbind, unname(chrom_list))

# 2. Separate by celltype
# no rownames (sample names) in files --> get them on a different way
pseudosamples <- read.csv(paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/Pseudosamples.csv"))
colnames(chrom_df) <- pseudosamples$Pseudosamples

# 3. Save matrix to file
write.csv(chrom_df, paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/Pseudobulk_GeneScoreMatrix", norm, "_AllChromosomes.csv"))