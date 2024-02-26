## Plot heatmap indicating how many cells of each samples
## are part of each cluster

# run with conda environment sc-atac

library(ArchR)
library(pheatmap)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")


# 2. Confusion matrix of clusters and samples

cM <- confusionMatrix(paste0(projPM$Clusters_r0.5), paste0(projPM$Sample))
cM

cM_harmony <- confusionMatrix(paste0(projPM$ClustHarmony_r0.5), paste0(projPM$Sample))
cM_harmony

# 3. Heatmap of confusion matrix

myPalette <- colorRampPalette(paletteContinuous("whiteBlue"), bias = 5)

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    #color = paletteContinuous("whiteBlue"), 
    #color = myPalette,
    border_color = "black"
)
p


cM_harmony <- cM_harmony / Matrix::rowSums(cM_harmony)
p_harmony <- pheatmap::pheatmap(
    mat = as.matrix(cM_harmony),
    #color = paletteContinuous("whiteBlue"),
    #color = myPalette,
    border_color = "black"
)
p_harmony

# 4. Save heatmap

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_pdf(
   p, 
   "ArchROutput/Plots/ClusterSample_Confusion.pdf",
   width=12,
   height=12
)


save_pheatmap_pdf(
   p_harmony,
   "ArchROutput/Plots/ClusterSample_Confusion_Harmony.pdf",
   width=12,
   height=12
)
