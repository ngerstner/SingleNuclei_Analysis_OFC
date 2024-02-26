## Plot UMAP visualization with cell type labels
## transferred from RNA data

# run with conda environment sc-atac

library(ArchR)
library(Seurat)
library(SeuratDisk)

# 1. Load ArchR project with ATACseq data
projPM <- loadArchRProject(path = "ArchROutput/")

# 1a. Load cell type labels created with In_RELN in reference
df_RELN <- readRDS(paste0(projPM@projectMetadata$outputDirectory, "/tables/07_1_Metadata_RELN.rds"))
projPM <- addCellColData(
  ArchRProj = projPM,
  data = df_RELN$predictedGroup_Un_RELN,
  name = "predictedGroup_Un_RELN",
  cells = rownames(df_RELN),
  force = FALSE
)


# 2. Plot cell type labels
p1 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "predictedGroup_Un"
)

p5 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "predictedGroup_Un_RELN"
)

p2 <- plotEmbedding(projPM,
	colorBy = "cellColData",
	name = "Clusters_r0.5"
)

p3 <- plotEmbedding(projPM,
        colorBy = "cellColData",
        name = "Clusters_r0.5"
)

p4 <- plotEmbedding(projPM,
        colorBy = "cellColData",
        name = "Clusters_r1"
)

# 3. Save plots as PDF
plotPDF(p1,p5,p2,p3,p4,
	name = "UMAP_IterativeLSI_cellTypes.pdf",
	ArchRProj = projPM,
	addDOC = FALSE,
	width = 8,
	height = 8)
