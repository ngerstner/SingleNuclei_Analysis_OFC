## Perform final filtering (one sample that clusters highly together and some clusters that have low quality -
## regarding doublet scores and number of fragments)  

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")

# 2. Filter out sample/clusters from ArchR object

filter_SU645 <- rownames(projPM@cellColData[projPM@cellColData$Sample != "SU645_PFC_ATAC",])
filter_clusters <- rownames(projPM@cellColData[projPM@cellColData$Clusters_r1 %ni% c("C8", "C19", "C20", "C21", "C28", "C29"),])

filter_intersect <- intersect(filter_SU645, filter_clusters)

# subset actual ArchRObject
projPM1 <- subsetArchRProject(projPM,
	cells = filter_intersect,
	outputDirectory = "ArchRSubset_FinalFiltering",
	force = TRUE,
	dropCells = FALSE)
