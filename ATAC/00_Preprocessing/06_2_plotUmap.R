## Plot UMAP representation with different colourings

# run with conda environment sc-atac

library(ArchR)
library(dplyr)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")

# add X6.Batch and Batch to cellColData
b <- inner_join(as.data.frame(projPM@cellColData), as.data.frame(projPM@sampleColData[c("ID", "Main.Batch")]), 
		by = c("Sample" = "ID"))
projPM <- addCellColData(projPM,
	       data = b$Main.Batch,
	       name = "Main.Batch",
	       cells = row.names(projPM@cellColData))
# save Batch as factor
projPM@cellColData$X6.Batch <- as.factor(projPM@cellColData$X6.Batch)
projPM@cellColData$Main.Batch <- as.factor(projPM@cellColData$Main.Batch)

# 2. Plot UMAP with different colourings

p1 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "Sample", embedding = "UMAP") 

p2 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "Main.Batch", embedding = "UMAP")

p3 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "X6.Batch", embedding = "UMAP")

p4 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "Clusters_r0.5", embedding = "UMAP")

p5 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "Clusters_r0.8", embedding = "UMAP")

p6 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "Clusters_r1", embedding = "UMAP")

p7 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")

p8 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

p9 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")

p10 <- plotEmbedding(ArchRProj = projPM, colorBy = "cellColData", name = "mito_perc", embedding = "UMAP")
# 3. Save plots in PDF

plotPDF(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
        name = "UMAPs_IterativeLSI.pdf",
        ArchRProj = projPM,
        addDOC = FALSE,
        width = 12,
        height = 15
)



# 2a. Plot UMAP with different colourings with smaller point size (somehow does not work)

p1 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "Main.Batch", embedding = "UMAP")

p3 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "X6.Batch", embedding = "UMAP")

p4 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "Clusters_r0.5", embedding = "UMAP")

p5 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "Clusters_r0.8", embedding = "UMAP")

p6 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "Clusters_r1", embedding = "UMAP")

p7 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")

p8 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

p9 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")

p10 <- plotEmbedding(ArchRProj = projPM, size = 0.0001, plotAs = "points", colorBy = "cellColData", name = "mito_perc", embedding = "UMAP")

# 3a. Save plots in PDF

plotPDF(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
        name = "UMAPs_IterativeLSI_pointSize0.0001.pdf",
        ArchRProj = projPM,
        addDOC = FALSE,
        width = 12,
        height = 15
)

