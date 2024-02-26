## TF motif enrichment analysis within the peaks nearby 
## previously identified DE/DA genes between PRS extreme groups

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(awtools)

source("07_1a_Functions_MotifEnrichment.R")

# suffix for alternative peak calling version
#alt_peak <- "_alt1"
alt_peak <- ""

# folder of ArchR project
# archr_folder <- "ArchR_Peaks"
archr_folder <- "../ArchRSubset_FinalFiltering"
basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# PRS trait
prs_trait <- "SCZ2022"
#prs_trait <- "crossDisorder2019"
#prs_trait <- "BIP2021"  # 0 DE DA hits
#prs_trait <- "MDD"
#prs_trait <- "height"   # 0 DE DA hits

# 1. Read data

projPM <- loadArchRProject(path = paste0(archr_folder, alt_peak ,"/"))

# 2. Add motif annotations to ArchRProject
# --> returns binary value for each peak indicating if a motif is observed within peak region

## Motifs from the chromVar "cisbp" dataset
# projPM <- addMotifAnnotations(ArchRProj = projPM, motifSet = "cisbp", annoName = "Motif")

# Motifs from the JASPAR2020 database
projPM <- addMotifAnnotations(ArchRProj = projPM, motifSet = "JASPAR2020", annoName = "Motif_JASPAR2020")

# 3. Analyze matches
#mo_matches <- getMatches(projPM, "Motif")
#str(mo_matches)

mo_matches <- getMatches(projPM, "Motif_JASPAR2020")
str(mo_matches)


# 4 For diff. acc. genes between cases and controls in each cell type

# 4.1 Read gene score results
p <- paste0("GeneScore_WaldTest_OutlierRemoved_Normalized_.*_", prs_trait, ".csv")
files <- list.files(path = paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/"),
                      pattern = p,
                      full.names = FALSE)
res_list <- list()
for (f in files){
	# read in results of differential gene score matrix analysis
	res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/", f))
    ct <- sub(paste0(".*_Normalized_(.*)_",prs_trait,"\\.csv"),"\\1",f)
	res_list[[ct]] <- res
}
res_df <- bind_rows(res_list, .id="celltype")

# 4.2 Read DE results
p <- paste0("sig0.1_Wald_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
              prs_trait, ".csv")
files <- list.files(path = paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/PRSext"),
                      pattern = p,
                      full.names = FALSE)
# read in the diff genes of different cell types
de_ct <- list()
for (f in files){
    de <- read.csv(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/PRSext/",f))
    ct <- sub("(\\w*)_~.*","\\1",f)
    de_ct[[ct]] <- de
}
de_ct <- de_ct[sapply(de_ct, 
            function(x) dim(x)[1]) > 0] # remove empty dataframes from list
celltypes_de <- names(de_ct)
df <- bind_rows(de_ct, .id="celltype")
df$celltype <- factor(df$celltype, levels = celltypes_de)
# subset to the significant hits
df$padj[is.na(df$padj)] <- 1
df_sig <- df[df$padj <= 0.1,] # TODO: --> FIX NAs here
head(df_sig)

# 4.3 FDR correction

# on merged object with DE hits
df_merge_de <- df_sig %>%
	inner_join(res_df, by=c("celltype", "X"="gene_name"))
de_df <- fdr_perCelltype(df_merge_de, "pvalue.y") 	
de_da_df <- de_df %>%
    filter(sig)

# for all genes
all_df <- fdr_perCelltype(res_df, "pvalue")
all_df_de <- all_df %>%
    filter(sig)
de_da_df <- all_df_de


# 4.4 Perform motif enrichment per celltype
enr_list <- list()
for (ct in unique(de_da_df$celltype)){

    # subset DE genes to respective celltype
    de_da_ct <- de_da_df %>%
        filter(celltype == ct)

    # get GRanges object for all peaks tested for motifs
    r <- rowRanges(mo_matches)
    # subset GRanges object to the peaks in promoter region of DE&DA genes
    matches_diff <- r %>%
        subset(nearestGene %in% de_da_ct$X)

    # perform enrichment analysis on respective peaks
    enr_df <- peakAnnoEnrichmentAdapted(
        seMarker = matches_diff,
        ArchRProj = projPM,
        peakAnnotation = "Motif_JASPAR2020",
        matches = mo_matches
    )

    enr_list[[ct]] <- enr_df[,c("feature", "mlog10Padj")]

}

# combine motif enrichment results from all celltypes in one dataframe
enr_all <- enr_list %>% purrr::reduce(full_join, by = "feature")
colnames(enr_all) <- c("feature", unique(de_da_df$celltype))
head(enr_all)

# filter to those motifs enriched in at least one cell type
enr_all_filt <- enr_all %>%
    filter_at(vars(2:ncol(enr_all)), any_vars(. >= -log10(0.05)))

# get TF class and family for each TF
tf_class <- sapply(enr_all_filt$feature, function(x){
    projPM@peakAnnotation$Motif_JASPAR2020$motifs[[x]]@matrixClass
})

tf_family <- sapply(enr_all_filt$feature, function(x){
    projPM@peakAnnotation$Motif_JASPAR2020$motifs[[x]]@tags$family
})


# 5. Plot enrichment results

# heatmap
mat <- as.matrix(enr_all_filt[,2:ncol(enr_all_filt)])
rownames(mat) <- enr_all_filt$feature
if(length(unique(tf_class)) <= 2){
    col_pal <- brewer.pal(n = 3, name="Set1")[1:length(unique(tf_class))]
} else {
    col_pal <- bpalette[1:length(unique(tf_class))]
}
ha <- rowAnnotation(TF_class = tf_class,
    col = list(TF_class = setNames(col_pal, unique(tf_class)))
    )
ht <- Heatmap(mat, 
    name = "-log10(P.adj)", 
    col = colorRamp2(breaks = c(0, 2, 4),
                    colors = brewer.pal(n=3, name="OrRd")),
    right_annotation = ha)

pdf(paste0(DEA_dir, "plots/07_1_MotifEnrichment/07_01_DEandDA_MotifEnrichment_PRSext_", prs_trait, ".pdf"),
    height=12, width=8)
pdf(paste0(DEA_dir, "plots/07_1_MotifEnrichment/07_01_DA_MotifEnrichment_PRSext_", prs_trait, ".pdf"),
    height=12, width=8)
draw(ht)#, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()

# r <- rowRanges(mo_matches)
# ino80e <- r %>%
#     subset(nearestGene == "INO80E")

# mo_matches_de <- mo_matches %>% subset(nearestGene == "INO80E")
# mo_matches_de@assays@data$matches

# n_matches <- apply(mo_matches_de@assays@data$matches, 2, sum)
# n_matches_min4 <- n_matches[n_matches>=4]

# Jaspar:
# KLF5_57 ZNF148_548   KLF4_587 
#       4          4          5 



