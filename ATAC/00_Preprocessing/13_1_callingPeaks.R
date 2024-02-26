## Calling peaks with MACS2 with default parameter

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)

# 1. Read data

projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")

# 2. Find path to MACS2

pathToMacs2 <- findMacs2()

# 3. Add peaks

projPM <- addReproduciblePeakSet(
    ArchRProj = projPM, 
    groupBy = "celltypes_refined", 
    pathToMacs2 = pathToMacs2
)

# 4. Write ArchR project with peaks to file
saveArchRProject(ArchRProj = projPM, outputDirectory = "ArchRSubset_FinalFiltering", load = FALSE)
