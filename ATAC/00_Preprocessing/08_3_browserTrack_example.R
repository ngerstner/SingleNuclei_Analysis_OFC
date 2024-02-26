## Browser Tracks for selected genes
## which were indicated in differential analysis

# run with conda environment sc-atac

library(ArchR)

# 1. Load ArchR project with ATACseq data
projPM <- loadArchRProject(path = "ArchROutput/")


# 2. Plot browser track for example gene

p <- plotBrowserTrack(
    ArchRProj = projPM, 
    groupBy = "celltypes_refined",  # refinement of celltypes only done in script 09_3
    geneSymbol = c("FKBP5", "TCF4", "NR3C1", "HTR2C", "ZBTB16", "HTR3A", "FAM87B", "IGFBP5"), 
    upstream = 100000,
    downstream = 100000
)

# 3. Save plots as files

plotPDF(plotList = p, 
    name = "08_3_Plot-Tracks-Examples.pdf", 
    ArchRProj = projPM, 
    addDOC = FALSE, width = 5, height = 5)