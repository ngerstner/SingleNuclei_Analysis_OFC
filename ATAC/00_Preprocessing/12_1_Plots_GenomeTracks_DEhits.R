## Genome Tracks for example genes indicated in DE analysis
## with one track per cell type

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
df_10 <- bind_rows(de_ct_10, .id="celltype")
df_10$celltype <- factor(df_10$celltype, levels = celltypes_de)

# 2. Track Plotting

# with gene integration matrix
p <- plotBrowserTrack(
    ArchRProj = projPM, 
    groupBy = "celltypes_refined", 
    geneSymbol = df_10$X, 
    useMatrix = "GeneIntegrationMatrix",
    upstream = 100000,
    downstream = 100000
)

plotPDF(plotList = p, 
    name = "10_7_GenomeTracks/Plot-Tracks-DEhits-max10perCelltype-GeneIntegration.pdf", 
    ArchRProj = projPM, 
    addDOC = FALSE, width = 5, height = 5)

# with gene score matrix
p <- plotBrowserTrack(
    ArchRProj = projPM, 
    groupBy = "celltypes_refined", 
    geneSymbol = df_10$X, 
    useMatrix = "GeneScoreMatrix",
    upstream = 100000,
    downstream = 100000
)

plotPDF(plotList = p, 
    name = "10_7_GenomeTracks/Plot-Tracks-DEhits-max10perCelltype-GeneScore.pdf", 
    ArchRProj = projPM, 
    addDOC = FALSE, width = 5, height = 5)

