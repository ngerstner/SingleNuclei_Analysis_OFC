## Genome Tracks for example genes indicated in DE analysis
## with one track per disease status

# run with conda environment sc-atac

## 0. Load packages and setup
library(ArchR)
library(dplyr)
library(RColorBrewer)
library(rlang) # necessary for function as_string

# define pathes to data and folders
# basedir <- '~/Documents/PostmortemBrain/workspace/'
basedir <- '/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/'

# 1. Read data

# 1.1 Read ArchR object with final filtering
projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")

# 1.2 Read DE hits
p <- "sig0.1_Wald_pcFromCorrectedDataAllCellTypes.csv"
files <- list.files(path = paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2/tables/DESeq2_Wald/"),
		pattern = p,
		full.names = FALSE)
de_ct <- list()
for (f in files){
	de <- read.csv(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2/tables/DESeq2_Wald/", f))
	ct <- sub("(\\w*)_~.*", "\\1", f)
	de_ct[[ct]] <- de
}
celltypes_de <- names(de_ct)
de_ct <- de_ct[sapply(de_ct, function(x) dim(x)[1]) > 0] # remove celltypes with zero hits

# list to dataframe
df <- bind_rows(de_ct, .id="celltype")
df$celltype <- factor(df$celltype, levels = celltypes_de)

# max 10 genes per celltype
de_ct_10 <- lapply(de_ct, function(x) x[sample(nrow(x), min(nrow(x), 10)),])


# 3. Track plotting

for (ct in names(de_ct_10)){

    # subset ArchR object to cells from celltype
    idxCelltype <- BiocGenerics::which(projPM$celltypes_refined %in% ct)
    cellsCelltype <- projPM$cellNames[idxCelltype]
    projCelltype <- projPM[cellsCelltype, ]

    # add column for disease status to cellColData
    cellMeta <- as.data.frame(projCelltype@cellColData)
    sampleMeta <- as.data.frame(projCelltype@sampleColData)
    sampleMeta$ID <- rownames(sampleMeta)
    cellMeta <- cellMeta %>%
        left_join(sampleMeta[c("ID", "Status")], by=c("Sample"="ID"))
    cellMeta$Status <- as.character(cellMeta$Status)
    projCelltype <- addCellColData(ArchRProj = projCelltype, data = cellMeta$Status,
        cells = projCelltype$cellNames, name = "Status")

    # with gene integration matrix
    p <- plotBrowserTrack(
        ArchRProj = projCelltype, 
        groupBy = "Status", 
        geneSymbol = de_ct_10[[ct]]$X, 
        useMatrix = "GeneIntegrationMatrix",
        upstream = 50000,
        downstream = 50000
    )

    plotPDF(plotList = p, 
        name = paste0("10_8_GenomeTracks_CaseControl/Plot-Tracks-DEhits-max10perCelltype-GeneIntegration_", ct, ".pdf"), 
        ArchRProj = projPM, 
        addDOC = FALSE, width = 5, height = 5)


    # with gene score matrix
    p <- plotBrowserTrack(
        ArchRProj = projCelltype, 
        groupBy = "Status", 
        geneSymbol = de_ct_10[[ct]]$X, 
        useMatrix = "GeneScoreMatrix",
        upstream = 50000,
        downstream = 50000
    )

    plotPDF(plotList = p, 
        name = paste0("10_8_GenomeTracks_CaseControl/Plot-Tracks-DEhits-max10perCelltype-GeneScore_", ct, ".pdf"), 
        ArchRProj = projPM, 
        addDOC = FALSE, width = 5, height = 5)



}


# p <- plotBrowserTrack(
#         ArchRProj = projCelltype, 
#         groupBy = "Status", 
#         geneSymbol = "HES4", 
#         useMatrix = "GeneScoreMatrix",
#         upstream = 100000,
#         downstream = 100000
#     )

#     plotPDF(plotList = p, 
#         name = paste0("10_8_GenomeTracks_CaseControl/Plot-Tracks-Example-Exc_L4-6_1_HES4-GeneScore.pdf"), 
#         ArchRProj = projPM, 
#         addDOC = FALSE, width = 5, height = 5)