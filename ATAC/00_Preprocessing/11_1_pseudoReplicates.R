## Generate pseudobulk replicates with ArchR method

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)

# 1. Read data

projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")


# 2. Generate group coverages/pseudobulk samples

#projPM <- addGroupCoverages(
#  ArchRProj = projPM,
#  groupBy = "celltypes_clusters_r1",
#  useLabels = TRUE,
#  minCells = 40,
#  maxCells = 500,
#  maxFragments = 25 * 10^6,
#  minReplicates = 2,
#  maxReplicates = 5,
#  sampleRatio = 0.8,
#  kmerLength = 6,
#  threads = getArchRThreads(),
#  returnGroups = FALSE,
#  parallelParam = NULL,
#  force = FALSE,
#  verbose = TRUE,
#  logFile = createLogFile("addGroupCoverages")
#)
projPM <- addGroupCoverages(
  ArchRProj = projPM,
  groupBy = "celltypes_refined",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 2000,
  maxFragments = 25 * 10^6,
  minReplicates = 90,
  maxReplicates = 90,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)


# 3. Save ArchR project for later usage
saveArchRProject(projPM)

