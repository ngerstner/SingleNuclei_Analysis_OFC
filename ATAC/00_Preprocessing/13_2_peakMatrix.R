## Add peak matrix to ArchR project and also save as rds-File

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

projPM <- loadArchRProject(path = paste0(archr_folder, alt_peak ,"/"))

# 2. Add peak matrix
projPM <- addPeakMatrix(projPM)
getAvailableMatrices(projPM)
mat_peaks <- getMatrixFromProject(projPM, useMatrix="PeakMatrix")
#mat_peaks@assays@data$PeakMatrix
saveRDS(mat_peaks, file = paste0(archr_folder, alt_peak, "/PeakMatrices/PeakMatrix", alt_peak, ".rds"))

# peak set
peak_set <- getPeakSet(projPM)

# 4. Write ArchR project with peaks to file
saveArchRProject(ArchRProj = projPM, outputDirectory = paste0(archr_folder, alt_peak), load = FALSE)