## Calculate correlation between gene scores and gene expression
## on normalized matrices (corrected for number of cells per pseudosample)
## across pseudosample (not averaged)

# run with conda environment sc-atac

library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(stringr)
library(org.Hs.eg.db)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# gene score matrix normalized?
norm <- "_Normalized"
# norm <- ""

celltypes <- c("Astro_FB", "Astro_PP", "Endothelial", "Exc_L2-3", "Exc_L3-5",
		"Exc_L4-6_1", "Exc_L4-6_2", "In_PVALB_Ba", "In_PVALB_Ch", 
		"In_SST", "In_VIP", "Microglia", "Oligodendrocyte", "OPC")


# 1. Read DE hits
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


# 2. Read gene expression pseudobulk data

# Read count 
counts <- readRDS(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2/data/counts.rds"))
rownames(counts) <- gsub("_RNA", "", rownames(counts))
counts[1:5,1:5]

# normalize count data
counts_norm <- t(counts)
counts_filt <- counts_norm[rowSums(counts_norm >= 5) >= 0.75*ncol(counts_norm), ]
vsd <- vst(as.matrix(counts_filt))
vsd <- t(vsd)
vsd[1:5,1:5]

# 2a. Split count data to list with entry per celltype
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(vsd), 
                                    pattern = "__",  
                                    n = 2), 
                 `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
vsd_split <- split.data.frame(vsd, 
                                 factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 str_split(rownames(u), '__', simplify = TRUE)[,2]))
# Explore the different components of list
str(vsd_split)


# 3. Read pseudo gene score matrix

geneScoreMat <- read.csv(paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/Pseudobulk_GeneScoreMatrix", norm, "_AllChromosomes.csv"),
					header = TRUE, row.names=1)
dim(geneScoreMat)

# correct sample IDs
colnames(geneScoreMat) <- gsub("(_ATAC\\.)", "#", colnames(geneScoreMat))
colnames(geneScoreMat) <- gsub("\\.", "-", colnames(geneScoreMat))
tmp <- str_split(colnames(geneScoreMat), pattern = "#", n = 2, simplify=TRUE)
colnames(geneScoreMat) <- paste0(tmp[,2], "__", tmp[,1])
geneScoreMat[1:5,1:5]

# filter genes in gene score matrix to make sure they are expressed in most of the samples
geneScoreMat_filt <- geneScoreMat[rowSums(geneScoreMat > 0) >= 0.75*ncol(geneScoreMat), ]
dim(geneScoreMat_filt)
geneScoreMat_filt[1:5,1:5]

# 3a. Split count data to list with entry per celltype
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(geneScoreMat_filt), 
                                    pattern = "__",  
                                    n = 2), 
                 `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
genescore_split <- split.data.frame(geneScoreMat_filt, 
                                 factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 str_split(rownames(u), '__', simplify = TRUE)[,2]))
# Explore the different components of list
str(genescore_split)


# 4. Get overall correlation across celltypes

counts <- vsd
geneScoreMat <- t(geneScoreMat_filt)

# reorder gene score matrix according to ordering of gene expression matrix
# and subset count matrix for common pseudosamples
indices <- intersect(rownames(counts), rownames(geneScoreMat))
geneScoreMat <- geneScoreMat[indices,]
counts <- counts[indices,]

# subset both matrices to common set of genes
colnames(geneScoreMat) <- mapIds(org.Hs.eg.db,
                     keys=colnames(geneScoreMat),
                     column="ENSEMBL",
                     keytype="SYMBOL",
                     multiVals="first")
indices_col <- intersect(colnames(counts), colnames(geneScoreMat))
geneScoreMat <- geneScoreMat[,indices_col]
counts <- counts[,indices_col]

counts[1:5,1:5]
geneScoreMat[1:5,1:5]

# get overall correlation
rowCor_mat <- sapply(1:ncol(counts), function(i) cor(counts[,i], geneScoreMat[,i], method = "pearson"))

# variance
var_exp <- sapply(1:ncol(counts), function(i) var(counts[,i]))
var_score <- sapply(1:ncol(counts), function(i) var(geneScoreMat[,i]))
# mean 
mean_exp <- sapply(1:ncol(counts), function(i) mean(counts[,i]))
mean_score <- sapply(1:ncol(counts), function(i) mean(geneScoreMat[,i]))

df_meanVar <- data.frame(gene_id = rep(colnames(counts), 2),
						modality = c(rep("gene expression", ncol(counts)),
									rep("gene score", ncol(counts))),
						variance = c(var_exp, var_score), 
						mean = c(mean_exp, mean_score))

# 4a. Plot correlation values
plot_df <- data.frame(gene_id = colnames(counts),
					corr = rowCor_mat)
plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
# without color scheme
ggplot(plot_df, aes(x=corr))+
	geom_histogram() +
	xlab("Correlation between gene expression and gene scores") +
	xlim(c(-1, 1)) +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationGeneScoresAndGeneExp_Normalized_Pseudobulk_Histogram.pdf"),
	width= 8,
	height=6)


# plot mean-variance scatterplot
ggplot(df_meanVar, aes(x=mean, y=variance, color=modality))+
	geom_point(alpha=0.1) + 
	xlab("Mean") +
	ylab("Variance") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationGeneScoresAndGeneExp_Normalized_Pseudobulk_MeanVariance.pdf"),
	width= 8,
	height=6)