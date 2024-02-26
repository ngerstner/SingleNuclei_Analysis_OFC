## Perform clustering with different resolutions (0.5, 0.8, 1)
## on Iterative LSI embedding

# run with conda environment sc-atac

library(ArchR)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")


# 2. Perform clustering

# on IterativeLSI embedding
projPM <- addClusters(
  input = projPM,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_r0.5",
  resolution = 0.5,
  maxClusters = 40
)

projPM <- addClusters(
  input = projPM,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_r0.8",
  resolution = 0.8,
  maxClusters = 40
)

projPM <- addClusters(
  input = projPM,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_r1",
  resolution = 1,
  maxClusters = 40
)

# on Harmony embedding
projPM <- addClusters(
  input = projPM,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustHarmony_r0.5",
  resolution = 0.5,
  maxClusters = 40
)

projPM <- addClusters(
  input = projPM,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustHarmony_r0.8",
  resolution = 0.8,
  maxClusters = 40
)

projPM <- addClusters(
  input = projPM,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "ClustHarmony_r1",
  resolution = 1,
  maxClusters = 40
)


# 3. Save project

saveArchRProject(
  ArchRProj = projPM,
  outputDirectory = "ArchROutput/",
)



