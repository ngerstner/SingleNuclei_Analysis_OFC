## Assign cell types by assigning the cell type to each cluster
## that the majority of nuclei are labeled

# run with conda environment sc-atac

library(ArchR)

# 1. Load ArchR project with ATACseq data
projPM <- loadArchRProject(path = "ArchROutput/")

# 2. Overlap between clusters and cell type labels
cM <- confusionMatrix(projPM$Clusters_r1, projPM$predictedGroup)

# 3. Assign cell type label with max overlap to each cluster
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]

projPM$celltypes_clusters_r1 <- mapLabels(projPM$Clusters_r1, newLabels = labelNew, oldLabels = labelOld)

# 4. Plot cell type labels
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



# 5. Save plots as PDF
plotPDF(p1,p2,p3,
	name = "UMAPs_IterativeLSI_cellTypes_clusters_r1.pdf",
	ArchRProj = projPM,
	addDOC = FALSE,
	width = 8,
	height = 9)


# 6. Save ArchR project
saveArchRProject(ArchRProj = projPM)
