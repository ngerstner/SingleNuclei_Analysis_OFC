## Plot genome tracks with ArchR for DE and DA hits between PRS extreme groups

# run with conda environment sc-atac

## 0. Load packages and setup
library(ArchR)
library(dplyr)
library(RColorBrewer)
library(rlang) # necessary for function as_string

# define pathes to data and folders
# basedir <- '~/Documents/PostmortemBrain/workspace/'
basedir <- '/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/'
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")


# 1. Read data

# 1.1 Read ArchR object with final filtering
projPM <- loadArchRProject(path = "../ArchRSubset_FinalFiltering/")

# 1.2 Define hits to plot
hits <- list()
hits[["SCZ2022"]] <- list()
hits[["SCZ2022"]][["Exc_L2-3"]] <- c("INO80E", "HCN2")


# 3. Track plotting

for (trait in names(hits)){

    # read samples in PRS extreme groups
    prsext_samples <- read.csv(paste0(DEA_dir, "tables/PRSext/samples_propensityScore/",
								trait, "_samples_propensityScore.csv")) %>%
					dplyr::rename(PRSext = ExtGroup) 

    for (ct in names(hits[[trait]])){

    # subset ArchR object to cells from celltype
    idxCelltype <- BiocGenerics::which(projPM$celltypes_refined %in% ct)
    cellsCelltype <- projPM$cellNames[idxCelltype]
    projCelltype <- projPM[cellsCelltype, ]

    # subset ArchR object to cells from samples of PRSext groups
    idxSamples <- BiocGenerics::which(projCelltype$Sample %in% prsext_samples$ID)
    cellsSamples <- projCelltype$cellNames[idxSamples]
    projSamples <- projCelltype[cellsSamples, ]

    # add column for PRSext group to cellColData
    cellMeta <- as.data.frame(projSamples@cellColData)
    cellMeta <- cellMeta %>%
        left_join(prsext_samples[c("ID", "PRSext")], by=c("Sample"="ID"))
    projSamples <- addCellColData(ArchRProj = projSamples, data = cellMeta$PRSext,
        cells = projSamples$cellNames, name = "PRSext")


    # with gene integration matrix
    p <- plotBrowserTrack(
        ArchRProj = projSamples, 
        groupBy = "PRSext", 
        geneSymbol = hits[[trait]][[ct]], 
        useMatrix = "GeneIntegrationMatrix",
        upstream = 50000,
        downstream = 50000
    )

    plotPDF(plotList = p, 
        name = paste0("DA_03_3f_GenomeTracks_PRSext/03_3f_GenomeTracks_DEandDAhits-GeneIntegration_", trait, "_", ct, ".pdf"), 
        ArchRProj = projPM, 
        addDOC = FALSE, width = 7, height = 5)


    # with gene score matrix
    p <- plotBrowserTrack(
        ArchRProj = projSamples, 
        groupBy = "PRSext", 
        geneSymbol = hits[[trait]][[ct]], 
        useMatrix = "GeneScoreMatrix",
        upstream = 50000,
        downstream = 50000
    )

    plotPDF(plotList = p, 
        name = paste0("DA_03_3f_GenomeTracks_PRSext/03_3f_GenomeTracks_DEandDAhits-GeneScore_", trait, "_", ct, ".pdf"), 
        ArchRProj = projPM, 
        addDOC = FALSE, width = 7, height = 5)

    }
}
