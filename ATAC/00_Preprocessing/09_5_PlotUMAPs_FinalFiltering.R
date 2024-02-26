## Plot UMAPs of final filtering (one sample that clusters highly together and some clusters that have low quality -
## regarding doublet scores and number of fragments)

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)

# 1. Read data

projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")

# 1a. Relabel In_LAMP5 to In_RELN
projPM@cellColData$celltypes_refined[projPM@cellColData$celltypes_refined == "In_LAMP5"] <- "In_RELN"

# 2. Plot cell type labels
p1 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "predictedGroup"
)

p2 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "Clusters_r1"
)

p3 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "celltypes_clusters_r1"
)

p4 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "celltypes_refined"
)




# 3. Save plots as PDF
plotPDF(p1,p2,p3,p4,
	name = "09_5_UMAPs_FinalFiltering_IterativeLSI_cellTypes_clusters_r1.pdf",
	ArchRProj = projPM,
	addDOC = FALSE,
	width = 8,
	height = 9)
