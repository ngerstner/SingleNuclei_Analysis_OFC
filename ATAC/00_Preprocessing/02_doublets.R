## Calculate doublet scores on snATAC-seq data

# INFO: use conda environment sc-atac to run this

# 0. load libraries
library(ArchR)


# 1. general settings
addArchRThreads(threads = 20)
addArchRGenome("hg38")


# 2. generate list with all arrow files
arrowfolder <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/ArrowFiles"

ArrowFiles <- list.files(path = arrowfolder, pattern = "\\.arrow", full.names=TRUE)


# 3. calculate doublet scores
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
