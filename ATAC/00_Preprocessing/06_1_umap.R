## Generate UMAP embedding based on IterativeLSI embedding
## without and with batch correction

# run with conda environment sc-atac

library(ArchR)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")


# 2. Generate UMAP embedding 

projPM <- addUMAP(
    ArchRProj = projPM, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

projPM <- addUMAP(
    ArchRProj = projPM,
    reducedDims = "Harmony",
    name = "UMAP_Harmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

# 3. Save project

saveArchRProject(
  ArchRProj = projPM,
  outputDirectory = "ArchROutput/",
)



