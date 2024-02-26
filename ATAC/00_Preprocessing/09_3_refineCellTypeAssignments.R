## Manually refine cell type assignments
## after having checked known marker genes in the UMAPs
## with In_RELN celltype

# run with conda environment sc-atac

library(ArchR)
library(Seurat)
library(SeuratDisk)

# 1. Read data
projPM <- loadArchRProject(path = "ArchROutput/")


# 2. Manual refinement of celltypes
projPM@cellColData$celltypes_refined <- as.vector(projPM@cellColData$celltypes_clusters_r1)
projPM@cellColData$celltypes_refined[projPM@cellColData$Clusters_r1 == "C17"] <- "In_RELN"


# 3. UMAP of refined cell types
p1 <- plotEmbedding(ArchRProj = projPM, 
	colorBy = "cellColData", 
	name = "celltypes_refined", 
	embedding = "UMAP")

plotPDF(p1, name = "09_3_Plot-UMAP-CelltypesRefined_RELN.pdf", ArchRProj = projPM, addDOC = FALSE, width = 5, height = 5)

# 4. Save ArchRProject
#saveArchRProject(projPM)
