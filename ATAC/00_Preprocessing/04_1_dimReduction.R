## Dimensionality reduction with iterative LSI
## and batch correction with harmony

# run with conda environment sc-atac

library(ArchR)
library(dplyr)

addArchRThreads(20)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")

# 2. Iterative LSI

projPM <- addIterativeLSI(
    ArchRProj = projPM,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 3, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(2), 
        sampleCells = 10000,
	max.clusters = 6, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

# 3. Add batch variable to cellColData

batch_var <- as.data.frame(projPM@sampleColData@listData)[c("ID","X6.Batch")]

cell_sample <- as.data.frame(projPM@cellColData@listData)["Sample"]

batch_joined <- inner_join(cell_sample, batch_var, by = c("Sample" = "ID"))

projPM@cellColData[,"X6.Batch"] <- batch_joined$X6.Batch

# 4. Run Harmony batch correction

projPM <- addHarmony(
  ArchRProj = projPM,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "X6.Batch"
)
## --> 4 iterations, warning: Quick-TRANSfer stage steps exceeded maximum (= 25889550)


# 5. Save project

saveArchRProject(
  ArchRProj = projPM,
  outputDirectory = "ArchROutput/",
)
