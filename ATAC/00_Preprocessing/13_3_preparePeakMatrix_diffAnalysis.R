## Prepare peak matrix to be used in differential analysis and 
## data integration, e.g. save sample metadata, save as mtx-File
## and generate pseudocounts per donor/cell type combination

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)

# suffix for alternative peak calling version
#alt_peak <- "_alt1"
alt_peak <- ""

# folder of ArchR project
# archr_folder <- "ArchR_Peaks"
archr_folder <- "ArchRSubset_FinalFiltering"

# 1. Read data

projPM <- loadArchRProject(path = paste0(archr_folder, alt_peak, "/"))

# save sampleColData as metadata object
saveRDS(projPM@sampleColData, paste0(archr_folder, alt_peak, "/PeakMatrices/SampleMetadata.rds"))

# 2. Read peak matrix
#mat_peaks@assays@data$PeakMatrix
#mat_peaks <- getMatrixFromProject(projPM, useMatrix="PeakMatrix")
#saveRDS(mat_peaks, paste0(archr_folder, alt_peak, "/PeakMatrices/PeakMatrix", alt_peak, ".rds"))
mat_peaks <- readRDS(paste0(archr_folder, alt_peak, "/PeakMatrices/PeakMatrix", alt_peak, ".rds"))

# save matrix in format that can be read by python
writeMM(mat_peaks@assays@data$PeakMatrix, file=paste0(archr_folder, alt_peak, "/PeakMatrices/PeakMatrix", alt_peak, ".mtx"))

# peak set
peak_set <- getPeakSet(projPM)
write.csv(peak_set, paste0(archr_folder, alt_peak, "/PeakMatrices/PeakSet", alt_peak, ".csv"),
	quote=FALSE)

# 3. Generate pseudocounts per celltype/sample combination

# add column with pseudobulk sample per cell
projPM@cellColData$pseudo_group <- paste0(projPM@cellColData$Sample, "#", projPM@cellColData$celltypes_clusters_r1)

# add pseudobulk id to peak matrix as colnames
mat_pseudo_prep <- mat_peaks@assays@data$PeakMatrix
colnames(mat_pseudo_prep) <- projPM@cellColData$pseudo_group

# sum up counts per pseudobulk sample
mat_pseudo_prep <- t(mat_pseudo_prep)
mat_pseudo <- rowsum(mat_pseudo_prep, row.names(mat_pseudo_prep))

# save pseudobulk matrix as rds file
saveRDS(mat_pseudo, paste0(archr_folder, alt_peak, "/PeakMatrices/Pseudobulk_PeakMatrix", alt_peak, ".rds"))

# save it also as csv to be able to read it with python
write.csv(mat_pseudo, paste0(archr_folder, alt_peak, "/PeakMatrices/Pseudobulk_PeakMatrix", alt_peak, ".csv"))


