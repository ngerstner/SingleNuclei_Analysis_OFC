## Calculate correlation between average gene scores and gene expression
## on normalized matrices (corrected for number of cells per pseudosample)

# run with conda environment DESeq2

library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(stringr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# gene score matrix normalized?
norm <- "_Normalized"
# norm <- ""

celltypes <- c("Astro_FB", "Astro_PP", "Endothelial", "Exc_L2-3", "Exc_L3-5",
		"Exc_L4-6_1", "Exc_L4-6_2", "In_PVALB_Ba", "In_PVALB_Ch", 
		"In_RELN", "In_SST", "In_VIP", "Microglia", "Oligodendrocyte", "OPC")


# 1. Read DE hits
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


# 2. Read gene expression pseudobulk data

# Read count 
counts <- readRDS(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/data/counts.rds"))
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
dim(geneScoreMat) #23963  1255

# correct sample IDs
colnames(geneScoreMat) <- gsub("(_ATAC\\.)", "#", colnames(geneScoreMat))
colnames(geneScoreMat) <- gsub("\\.", "-", colnames(geneScoreMat))
tmp <- str_split(colnames(geneScoreMat), pattern = "#", n = 2, simplify=TRUE)
colnames(geneScoreMat) <- paste0(tmp[,2], "__", tmp[,1])
geneScoreMat[1:5,1:5]

# filter genes in gene score matrix to make sure they are expressed in most of the samples
# this is not counts data, but gene scores
geneScoreMat_filt <- geneScoreMat[rowSums(geneScoreMat > 0.1) >= 0.5*ncol(geneScoreMat), ]
dim(geneScoreMat_filt) #20242  1255
geneScoreMat_filt[1:5,1:5]
geneScoreMat_filt <- t(geneScoreMat_filt)

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
geneScoreMat <- geneScoreMat_filt

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

## A. Average gene expression/gene score

# get average expression/score per gene
exp_mean <- colMeans(counts)
score_mean <- colMeans(geneScoreMat)

# get overall correlation
avg_cor <- cor(exp_mean, score_mean)
print("All celltypes")
print(avg_cor)

# 4a. Plot correlation values
plot_df <- data.frame(gene_id = names(exp_mean),
					expression = exp_mean,
					score = score_mean)
plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
# all genes, colored according to DE hit (TRUE/FALSE)
ggplot(plot_df, aes(x=expression, y=score))+
	geom_point(aes(color=de_hit), alpha=0.5) +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	stat_cor(method="pearson") +
	xlab("Gene expression (RNA-seq)") +
	ylab("Gene score (ATAC-seq)") +
	ggtitle("Correlation between average gene expression and gene scores") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationAvgGeneScoresAndGeneExp_Normalized_Pseudobulk.pdf"),
	width= 8,
	height=6)

## B. Variance of gene expression/gene score

# get variance of expression/score per gene
exp_var <- colVars(counts)
score_var <- colVars(geneScoreMat)

# get overall correlation between variances
avg_cor <- cor(exp_var, score_var)
print("All celltypes")
print(avg_cor)

# 4a. Plot correlation values
plot_df <- data.frame(gene_id = names(exp_var),
					expression = exp_var,
					score = score_var)
plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
# all genes
ggplot(plot_df, aes(x=expression, y=score))+
	#geom_point(aes(color=de_hit), alpha=0.5) +
	geom_point(alpha=0.5) +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	stat_cor(method="pearson") +
	xlab("Variance of gene expression (RNA-seq)") +
	ylab("Variance of gene score (ATAC-seq)") +
	ggtitle("Correlation between variance of gene expression and gene scores") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationVarGeneScoresAndGeneExp_Normalized_Pseudobulk.pdf"),
	width= 8,
	height=6)
	

# 4b. Get correlation per celltype
for (ct in celltypes){

	if (ct %in% names(genescore_split)){

		counts_ct <- t(vsd_split[[ct]])
		genescore_ct <- t(genescore_split[[ct]])

		# subset both matrices to common set of genes
		colnames(genescore_ct) <- mapIds(org.Hs.eg.db,
							keys=colnames(genescore_ct),
							column="ENSEMBL",
							keytype="SYMBOL",
							multiVals="first")
		indices_col <- intersect(colnames(counts_ct), colnames(genescore_ct))
		genescore_ct <- genescore_ct[,indices_col]
		counts_ct <- counts_ct[,indices_col]

		counts_ct[1:5,1:5]
		genescore_ct[1:5,1:5]

		## A. Average gene expression/gene score

		# get average expression/score per gene
		exp_mean <- colMeans(counts_ct)
		score_mean <- colMeans(genescore_ct)

		# get overall correlation
		avg_cor <- cor(exp_mean, score_mean)
		print(ct)
		print(avg_cor)
	
		# Plot correlation values
		plot_df <- data.frame(gene_id = names(exp_mean),
							expression = exp_mean,
							score = score_mean)
		plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
		
		# all genes, colored according to DE hit (TRUE/FALSE)
		ggplot(plot_df, aes(x=expression, y=score))+
			geom_point(aes(color=de_hit), alpha=0.5) +
			geom_smooth(method = "lm", color = "black", se = FALSE) +
			stat_cor(method="pearson") +
			xlab("Gene expression (RNA-seq)") +
			ylab("Gene score (ATAC-seq)") +
			ggtitle(paste0(ct, ": Correlation between average gene expression and gene scores")) +
			theme_light() +
			theme(
				axis.title.x = element_text(size=15),
				axis.title.y = element_text(size=15),
				axis.text.x = element_text(size=12),
				axis.text.y = element_text(size=12),
				legend.text = element_text(size=12)
			)
		ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationAvgGeneScoresAndGeneExp_Pseudobulk_",ct,".pdf"),
			width= 8,
			height=6)


		## B. Variance of gene expression/gene score

		# get variance of expression/score per gene
		exp_var <- colVars(counts_ct)
		score_var <- colVars(genescore_ct)

		# get overall correlation between variances
		avg_cor <- cor(exp_var, score_var)
		print(ct)
		print(avg_cor)

		# 4a. Plot correlation values
		plot_df <- data.frame(gene_id = names(exp_var),
							expression = exp_var,
							score = score_var)
		plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
		# all genes
		ggplot(plot_df, aes(x=expression, y=score))+
			#geom_point(aes(color=de_hit), alpha=0.5) +
			geom_point(alpha=0.5) +
			geom_smooth(method = "lm", color = "black", se = FALSE) +
			stat_cor(method="pearson") +
			xlab("Variance of gene expression (RNA-seq)") +
			ylab("Variance of gene score (ATAC-seq)") +
			ggtitle(paste0(ct,": Correlation between variance of gene expression and gene scores")) +
			theme_light() +
			theme(
				axis.title.x = element_text(size=15),
				axis.title.y = element_text(size=15),
				axis.text.x = element_text(size=12),
				axis.text.y = element_text(size=12),
				legend.text = element_text(size=12)
			)
		ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_CorrelationVarGeneScoresAndGeneExp_Normalized_Pseudobulk_",ct, ".pdf"),
			width= 8,
			height=6)
			}
}


# 5. Random correlations

# get average expression/score per gene
exp_mean <- colMeans(counts)
score_mean <- sample(colMeans(geneScoreMat))	# reorder randomly

# get overall correlation
avg_cor <- cor(exp_mean, score_mean)
print("All celltypes")
print(avg_cor)

# 5a. Plot correlation values
plot_df <- data.frame(gene_id = names(exp_mean),
					expression = exp_mean,
					score = score_mean)
plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
# all genes, colored according to DE hit (TRUE/FALSE)
ggplot(plot_df, aes(x=expression, y=score))+
	geom_point(aes(color=de_hit), alpha=0.5) +
	geom_smooth(method = "lm", color = "black", se = FALSE) +
	stat_cor(method="pearson") +
	xlab("Gene expression (RNA-seq)") +
	ylab("Gene score (ATAC-seq)") +
	ggtitle("Correlation between average gene expression and gene scores") +
	theme_light() +
	theme(
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15),
		axis.text.x = element_text(size=12),
		axis.text.y = element_text(size=12),
		legend.text = element_text(size=12)
	)
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_RandomCorrelationAvgGeneScoresAndGeneExp_Normalized_Pseudobulk.pdf"),
	width= 8,
	height=6)


# 5b. Get correlation per celltype
for (ct in celltypes){

	if (ct %in% names(genescore_split)){

		counts_ct <- t(vsd_split[[ct]])
		genescore_ct <- t(genescore_split[[ct]])

		# subset both matrices to common set of genes
		colnames(genescore_ct) <- mapIds(org.Hs.eg.db,
							keys=colnames(genescore_ct),
							column="ENSEMBL",
							keytype="SYMBOL",
							multiVals="first")
		indices_col <- intersect(colnames(counts_ct), colnames(genescore_ct))
		genescore_ct <- genescore_ct[,indices_col]
		counts_ct <- counts_ct[,indices_col]

		counts_ct[1:5,1:5]
		genescore_ct[1:5,1:5]

		# get average expression/score per gene
		exp_mean <- colMeans(counts_ct)
		score_mean <- sample(colMeans(genescore_ct)) # reorder randomly

		# get overall correlation
		avg_cor <- cor(exp_mean, score_mean)
		print(ct)
		print(avg_cor)
	
		# Plot correlation values
		plot_df <- data.frame(gene_id = names(exp_mean),
							expression = exp_mean,
							score = score_mean)
		plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
		
		# all genes, colored according to DE hit (TRUE/FALSE)
		ggplot(plot_df, aes(x=expression, y=score))+
			geom_point(aes(color=de_hit), alpha=0.5) +
			geom_smooth(method = "lm", color = "black", se = FALSE) +
			stat_cor(method="pearson") +
			xlab("Gene expression (RNA-seq)") +
			ylab("Gene score (ATAC-seq)") +
			ggtitle(paste0(ct, ": Correlation between average gene expression and gene scores")) +
			theme_light() +
			theme(
				axis.title.x = element_text(size=15),
				axis.title.y = element_text(size=15),
				axis.text.x = element_text(size=12),
				axis.text.y = element_text(size=12),
				legend.text = element_text(size=12)
			)
		ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationAvgGeneScoreAndGeneExp_Pseudobulk/03_3b_RandomCorrelationAvgGeneScoresAndGeneExp_Pseudobulk_",ct,".pdf"),
			width= 8,
			height=6)
	}
}

