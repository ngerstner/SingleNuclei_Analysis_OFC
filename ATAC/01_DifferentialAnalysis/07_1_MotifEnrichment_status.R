## TF motif enrichment analysis within the peaks nearby 
## previously identified DE/DA genes between cases and controls

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

source("07_1a_Functions_MotifEnrichment.R")

# suffix for alternative peak calling version
#alt_peak <- "_alt1"
alt_peak <- ""

# folder of ArchR project
# archr_folder <- "ArchR_Peaks"
archr_folder <- "../ArchRSubset_FinalFiltering"
basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

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
head(str(mo_matches))


# 4 For diff. acc. genes between cases and controls in each cell type

# 4.1 Read gene score results
res_list <- list()
for (ct in unique(projPM@cellColData$celltypes_refined)){
	p <- paste0("GeneScore_WaldTest_OutlierRemoved_Normalized_", ct, ".csv")
	# read in results of differential gene score matrix analysis
	res <- fread(paste0("tables/GeneScoreMatrix/", p))
	res_list[[ct]] <- res
}
res_df <- bind_rows(res_list, .id="celltype")

# 4.2 Read DE results
p <- "X6.Batch_Wald_pcFromCorrectedDataAllCellTypes.csv"
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

# # for all genes
# all_df <- fdr_perCelltype(res_df, "pvalue")
# all_df_de <- all_df %>%
#     filter(sig)
# de_da_df <- all_df_de


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
ha <- rowAnnotation(TF_class = tf_class,
    col = list(TF_class = setNames(brewer.pal(n = length(unique(tf_class)), name = "Set1"), unique(tf_class)))
    )
ht <- Heatmap(mat, 
    name = "-log10(P.adj)", 
    col = colorRamp2(breaks = c(0, 2, 4),
                    colors = brewer.pal(n=3, name="OrRd")),
    right_annotation = ha)
pdf(paste0(DEA_dir, "plots/07_1_MotifEnrichment/07_01_DEandDA_MotifEnrichment_status.pdf"))
pdf(paste0(DEA_dir, "plots/07_1_MotifEnrichment/07_01_DA_MotifEnrichment_status.pdf"))#
draw(ht)#, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()



# 6. Enrichment analysis for INO80E only
gene_interest <- "INO80E"
# gene_interest <- "HCN2"

r <- rowRanges(mo_matches)
ino80e <- r %>%
    subset(nearestGene == gene_interest)

mo_matches_de <- mo_matches %>% subset(nearestGene == gene_interest)
colnames(mo_matches_de@assays@data$matches)
mo_matches_de@assays

n_matches <- apply(mo_matches_de@assays@data$matches, 2, sum)
n_matches_min4 <- n_matches[n_matches>=4]

# Jaspar:
# KLF5_57 ZNF148_548   KLF4_587 
#       4          4          5 

# perform enrichment analysis on respective peaks
enr_df <- peakAnnoEnrichmentAdapted(
    seMarker = ino80e,
    ArchRProj = projPM,
    peakAnnotation = "Motif_JASPAR2020",
    matches = mo_matches
)

enr_df[,c("feature", "mlog10Padj")]
# filter to those motifs enriched in at least one cell type
enr_ino80e_filt <- enr_df[,c("feature", "mlog10Padj")] %>%
    filter(mlog10Padj >= -log10(0.05))

enr_ino80e_filt
write.csv(enr_ino80e_filt, 
        paste0(DEA_dir, "tables/MotifEnrichment/07_1_MotifEnrichment_", gene_interest, ".csv"))
#          feature mlog10Padj
# KLF4_587 KLF4_587   2.080057


# HCN2:

#                          feature mlog10Padj
# MAZ_410                   MAZ_410  11.096357
# ZNF148_548             ZNF148_548  10.267757
# ZNF740_525             ZNF740_525   7.689057
# EGR1_563                 EGR1_563   7.625157
# KLF15_401               KLF15_401   6.955257
# KLF4_587                 KLF4_587   6.564757
# KLF11_400               KLF11_400   6.380857
# ZBTB14_545             ZBTB14_545   6.210657
# CTCFL_560               CTCFL_560   6.069157
# ZNF460_476             ZNF460_476   5.581957
# KLF5_57                   KLF5_57   5.420357
# TFAP2A.var.2_231 TFAP2A.var.2_231   5.121057
# TFAP2C.var.2_627 TFAP2C.var.2_627   4.942257
# SP9_449                   SP9_449   4.781157
# KLF9_588                 KLF9_588   4.583957
# ASCL1.var.2_526   ASCL1.var.2_526   4.504257
# ZIC5_468                 ZIC5_468   4.376557
# TFAP2B_232             TFAP2B_232   4.280557
# KLF16_169               KLF16_169   4.114957
# EGR3_162                 EGR3_162   3.987757
# TFAP2E_454             TFAP2E_454   3.971057
# TFDP1_318               TFDP1_318   3.885657
# PLAG1_28                 PLAG1_28   3.879457
# RREB1_11                 RREB1_11   3.325557
# TCF4_622                 TCF4_622   3.225157
# TFAP2A_626             TFAP2A_626   3.141757
# E2F6_561                 E2F6_561   3.033657
# SP3_522                   SP3_522   2.962257
# SP4_117                   SP4_117   2.789457
# ZIC4_172                 ZIC4_172   2.761557
# TFAP2C_235             TFAP2C_235   2.694457
# SP8_170                   SP8_170   2.654857
# EBF1_562                 EBF1_562   2.534657
# ZEB1_320                 ZEB1_320   2.518957
# TCF12.var.2_543   TCF12.var.2_543   2.449557
# EBF3_532                 EBF3_532   2.301857
# KLF14_168               KLF14_168   2.116057
# KLF17_402               KLF17_402   2.089657
# NRF1_40                   NRF1_40   1.995457
# ESR2_47                   ESR2_47   1.925057
# EGR4_163                 EGR4_163   1.888357
# EGR2_161                 EGR2_161   1.871657
# ZNF135_470             ZNF135_470   1.854157
# KLF10_399               KLF10_399   1.799957
# SP2_521                   SP2_521   1.795057
# SP1_520                   SP1_520   1.604557
# ZNF263_632             ZNF263_632   1.553157
# ZBTB7B_128             ZBTB7B_128   1.542657
# ASCL1_482               ASCL1_482   1.536557
# TFAP2B.var.2_233 TFAP2B.var.2_233   1.429457
# ZBTB7C_129             ZBTB7C_129   1.327757
# GLIS2_165               GLIS2_165   1.326757
# TCF3_621                 TCF3_621   1.324957
# TCFL5_523               TCFL5_523   1.315657