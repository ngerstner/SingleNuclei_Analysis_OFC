## Extract pseudobulk gene score matrices per chromosome
## and save as csv-Files, 20.06.2023

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)
library(aod)
library(data.table)

# get argument for chromosome
args = commandArgs(TRUE)
chrom = args[1]
#chrom = "chr1"
#chrom = "chr2"
#chrom = "chr8"

# 1. Read data

projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")

# 2. Get gene score matrix
gsm <- getMatrixFromProject(
  ArchRProj = projPM,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = c(chrom),
  threads = 1
)

# extract actual matrix
mat_genescore <- gsm@assays@data$GeneScoreMatrix
rownames(mat_genescore) <- rowData(gsm)$name
mat_genescore[1:5,1:5]

# 3. Write gene score matrix to file
saveRDS(mat_genescore, file=paste0(projPM@projectMetadata$outputDirectory, "/GeneScoreMatrices/",
    "GeneScoreMatrix_", chrom, ".rds"))