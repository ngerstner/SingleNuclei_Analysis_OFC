## Analysis of correlation of gene expression and nearby peaks
## in addition to the number of peaks nearby genes

# run with conda environment DESeq2

## 0. Load packages
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(Hmisc)
library(org.Hs.eg.db)
library(ggpubr)


## 1. Setup

# define pathes to data and folders
# basedir <- '~/Documents/PostmortemBrain/workspace/'
basedir <- '/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/'
atacdir <- paste0(basedir, 'scripts/ATAC/')
# tabledir <- paste0(atacdir, 'ArchRSubset_FinalFiltering/PeakMatrices/')
tabledir <- paste0(atacdir, 'ArchR_Peaks/PeakMatrices/')
DEA_dir <- paste0(atacdir, '02_DifferentialAnalysis/')

# QC/Filtering
filt <- "_newFiltering"
# filt <- ""

# celltype of interest
args = commandArgs(TRUE)
cluster = args[1]
# cluster <- "Exc_L2-3"
# cluster <- "Oligodendrocyte"
# cluster <- "OPC"
# cluster <- "Astro_FB"
# cluster <- "Exc_L3-5"
# cluster <- "Exc_L4-6_2"
# cluster <- "In_RELN"


### 1. Read data

# 1.1 Read genes with their genomic coordinates
# CAUTION: these coordinates are different than the ones I downloaded for the scGLUE processing
# not clear why, but both are hg38
gene_anno <- read.csv(paste0(DEA_dir, "tables/GeneAnnotations/gene_annotations.csv"))
gene_anno$symbol[is.na(gene_anno$symbol)] <- ""
gene_anno$ensembl <- mapIds(org.Hs.eg.db,
                                  keys = as.character(gene_anno$gene_id),
                                  column = "ENSEMBL",
                                  keytype = "ENTREZID")
gene_anno$ensembl[is.na(gene_anno$ensembl)] <- ""
TSS_anno <- read.csv(paste0(DEA_dir, "tables/GeneAnnotations/TSS_annotations.csv"))

# 1.2 Read metadata of peaks
peaks <- read.csv(paste0(tabledir, "PeakSet.csv")) %>%
  separate(GroupReplicate, c("celltype", "sample"), sep = "\\._\\.") 
peaks$celltype <- str_replace(peaks$celltype, "LAMP5", "RELN")
# assign a peak id as colname
peaks$ID <- paste0(peaks$seqnames,"#",peaks$start,"#",peaks$end,"#",peaks$celltype)
#peaks <- peaks[!is.na(peaks$nearestGene),]
# CAUTION: annotation of nearest TSS name (tx_name) is wrong
# each peak has a width of 501 basepairs

# add column with ensembl id of nearest gene
peaks$nearestGene_ensembl <- mapIds(org.Hs.eg.db,
                                  keys = peaks$nearestGene,
                                  column = "ENSEMBL",
                                  keytype = "SYMBOL")

# subset peaks to promoter peaks with "distToTSS" <= 1000
peaks_promoter <- peaks %>%
    filter(distToTSS <= 1000) %>%
    filter(peakType == "Promoter")


# 1.2 Read counts and sample metadata
# count data
peak_counts <- readRDS(paste0(tabledir, "Pseudobulk_PeakMatrix.rds"))
colnames(peak_counts) <- peaks$ID
rownames(peak_counts) <- gsub("_ATAC", "", rownames(peak_counts))
peak_counts[1:5,1:5]

# subset counts to pseudosamples from specified cell type
splitted_ids_peaks <- as.data.frame(str_split(rownames(peak_counts), '#', simplify = TRUE))
samples_cluster_peaks <- which(splitted_ids_peaks$V2 == cluster)
rownames(peak_counts) <- splitted_ids_peaks$V1
cluster_peak_counts <- peak_counts[samples_cluster_peaks,]
# dim(cluster_peak_counts) # zB 90 742624

# normalize count data
counts_peak_norm <- t(cluster_peak_counts)
counts_peak_filt <- counts_peak_norm[rowSums(counts_peak_norm >= 5) >= 0.5*ncol(counts_peak_norm), ]
peak_vsd <- vst(as.matrix(counts_peak_filt))
peak_vsd <- t(peak_vsd)
peak_vsd[1:5,1:5]

# 1.3 Read gene expression pseudobulk data

# Read count 
counts <- as.matrix(readRDS(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/data/counts.rds")))
rownames(counts) <- gsub("_RNA", "", rownames(counts))
counts[1:5,1:5]

# subset counts to pseudosamples from specified cell type
splitted_ids <- as.data.frame(str_split(rownames(counts), '__', simplify = TRUE))
samples_cluster <- which(splitted_ids$V1 == cluster)
rownames(counts) <- splitted_ids$V2
cluster_counts <- counts[samples_cluster,]
# dim(cluster_counts) # zB 87 26195

# normalize count data
counts_norm <- t(cluster_counts)
counts_filt <- counts_norm[rowSums(counts_norm >= 5) >= 0.75*ncol(counts_norm), ]
vsd <- vst(as.matrix(counts_filt))
vsd <- t(vsd)
vsd[1:5,1:5]

# 1.4 Subset RNA and ATAC data to same set of samples
sample_intersect <- intersect(rownames(peak_vsd), rownames(vsd))
peak_vsd <- peak_vsd[sample_intersect,]
vsd <- vsd[sample_intersect,]

# 1.5 Sample metadata

# metadata
metadata <- readRDS(paste0(tabledir, "SampleMetadata.rds")) %>%
  as.data.frame %>%
  mutate_at(c("Sex", "Status", "X6.Batch"), factor)
metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))
metadata$RIN <- with(metadata, impute(RIN, median))
head(metadata)
colnames(metadata)

# 1.6 Read DE hits
p <- "sig0.1_Wald_pcFromCorrectedDataAllCellTypes.csv"
files <- list.files(path = paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/"),
		pattern = p,
		full.names = FALSE)
de_ct <- list()
for (f in files){
	de <- read.csv(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/", f))
	ct <- sub("(\\w*)_~.*", "\\1", f)
	de_ct[[ct]] <- de
}
celltypes_de <- names(de_ct)
de_ct <- de_ct[sapply(de_ct, function(x) dim(x)[1]) > 0] # remove celltypes with zero hits

# list to dataframe
df <- bind_rows(de_ct, .id="celltype")
df$celltype <- factor(df$celltype, levels = celltypes_de)

# 2. Correlation between gene expression and nearby genes

corr_gene_peaks <- function(gene_id, rna_mat, atac_mat, peak_meta){

  # find peaks with respective gene as nearest gene
  anno_peaks <- gene_anno[gene_anno$ensembl == gene_id,] %>%
      left_join(peaks, by="seqnames") %>%
      filter(start.y >= (start.x-100000), end.y <= (end.x+100000)) %>%
      filter(ensembl %in% colnames(vsd)) %>%
      filter(ID %in% colnames(peak_vsd))

  # calculate correlation between gene expression and accessibility data of each peak
  res_tab <- sapply(anno_peaks$ID, function(x){
                        cor_gene <- cor.test(vsd[,gene_id], peak_vsd[,x])
                        return(c("corr" = cor_gene$estimate,
                                "pvalue" = cor_gene$p.value))
                        },
                        simplify=FALSE, USE.NAMES=TRUE)

  return(res_tab)
}

# apply function to each gene in expression data
corr_peaks_per_gene <- sapply(colnames(vsd), corr_gene_peaks, vsd, peak_vsd, peaks, simplify = FALSE, USE.NAMES=TRUE)
corr_peak_df <- bind_rows(unlist(corr_peaks_per_gene, recursive=FALSE), .id = "gene.peak") %>%
      separate(gene.peak, c("gene", "peak"), "\\.") %>%
      mutate(fdr = p.adjust(pvalue, method="fdr"))  # 5 peaks with fdr <= 0.05

# number of peaks tested within window per gene
n_peaks_per_gene <- unlist(sapply(corr_peaks_per_gene, function(x) length(x)))
df_n_peaks <- data.frame("gene" = colnames(vsd),
                        "number_peaks" = n_peaks_per_gene)
# peaks tested for de genes (10 de genes)
df_n_peaks_de <- df_n_peaks %>%
    filter(gene %in% de_ct[[ct]]$gene)
str(df_n_peaks_de)

# histogram of number of peaks per gene
ggplot(df_n_peaks, aes(x=number_peaks))+
	geom_histogram() +
	xlab("Number of peaks within +/- 100000 bp from gene body") +
	ylab("Frequency") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3e_CorrelationRNAandPeak/03_3e_NumberOfPeaksPerGene_", cluster, "_Histogram.pdf"),
	width= 8,
	height=6)

# number of significantly correlated peaks per gene
head(sapply(corr_peaks_per_gene, function(x) sapply(x, "[[", 2)))
n_peaks_sig <- unlist(sapply(corr_peaks_per_gene, function(x) length(which(sapply(x, "[[", 2) <= 0.05))))
df_n_sig_peaks <- data.frame("gene" = colnames(vsd),
                        "number_peaks" = n_peaks_sig)

# histogram of (nominal) sig. correlated number of peaks per gene
ggplot(df_n_sig_peaks, aes(x=number_peaks))+
	geom_histogram(bins = max(df_n_sig_peaks$number_peaks)+1) +
	xlab("Number of sig. correlated peaks within +/- 100000 bp from gene body") +
	ylab("Frequency") +
  scale_y_log10() +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3e_CorrelationRNAandPeak/03_3e_NumberOfSigCorrPeaksPerGene_", cluster, "_Histogram.pdf"),
	width= 8,
	height=6)


# number of (fdr corrected) sig peaks per gene when only testing DE genes
corr_peak_df_de <- bind_rows(unlist(corr_peaks_per_gene, recursive=FALSE), .id = "gene.peak") %>%
      separate(gene.peak, c("gene", "peak"), "\\.") %>%
      filter(gene %in% de_ct[[ct]]$gene) %>%    # 79 peaks for 10 DE genes
      mutate(fdr = p.adjust(pvalue, method="fdr"))
head(corr_peak_df_de)
n_fdr <- length(which(corr_peak_df_de$fdr <= 0.05))
n_nominal <- length(which(corr_peak_df_de$pvalue <= 0.05))

# histogram of correlations for peaks around DE hits
ggplot(corr_peak_df_de, aes(x=corr.cor))+
	geom_histogram() +
	xlab("Correlation between gene expression and chromatin accessibility (+/- 100000 bp from gene body)") +
  #scale_y_log10() +
  ggtitle(paste0("Correlations between DE hits and nearby peaks: ", n_nominal, " nominally significant, ", n_fdr, " FDR significant")) +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3e_CorrelationRNAandPeak/03_3e_CorrelationsDEHitsAndPeaks_", cluster, "_Histogram.pdf"),
	width= 8,
	height=6)


# 3. Plot expression levels for genes against number of peaks

# get mean expression level per gene across samples
mean_exp <- as.data.frame(apply(vsd, 2, mean))
mean_exp$ensembl <- rownames(mean_exp)
colnames(mean_exp) <- c("meanExp", "ensembl")

# join mean expression with number of peaks and corr. peaks
df_meanExp <- inner_join(mean_exp, df_n_peaks, by=c("ensembl"="gene"))
df_meanExp <- inner_join(df_meanExp, df_n_sig_peaks, by=c("ensembl" = "gene"))
head(df_meanExp)

ggplot(df_meanExp, aes(x=number_peaks.x, y=meanExp))+
	geom_point(color="darkblue") +
	xlab("Number of peaks within +/- 100000 bp from gene body") +
	ylab("Mean expression level") +
	geom_smooth(method = "lm", color = "darkred", se = FALSE) +
	stat_cor(method="pearson") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3e_CorrelationRNAandPeak/03_3e_MeanExpression_vs_NumberOfPeaks_", cluster, ".pdf"),
	width= 7,
	height=7)

ggplot(df_meanExp, aes(x=number_peaks.y, y=meanExp))+
	geom_point(color="darkblue") +
	xlab("Number of sig. correlated peaks within +/- 100000 bp from gene body") +
	ylab("Mean expression level") +
	geom_smooth(method = "lm", color = "darkred", se = FALSE) +
	stat_cor(method="pearson") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3e_CorrelationRNAandPeak/03_3e_MeanExpression_vs_NumberOfCorrPeaks_", cluster, ".pdf"),
	width= 7,
	height=7)


