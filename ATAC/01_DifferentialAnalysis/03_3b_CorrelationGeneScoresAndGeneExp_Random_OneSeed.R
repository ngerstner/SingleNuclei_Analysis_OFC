## Correlation between gene scores and gene expression
## compare to random correlation (just one seed in this version)

# run with conda environment sc-atac

library(dplyr)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(stringr)
library(org.Hs.eg.db)
library(ggplot2)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

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

# 2a. Split count data to list with entry per celltype
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(counts), 
                                    pattern = "__",  
                                    n = 2), 
                 `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
counts_split <- split.data.frame(counts, 
                                 factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 str_split(rownames(u), '__', simplify = TRUE)[,2]))
# Explore the different components of list
str(counts_split)

# Subset the counts to only the specific celltype
#counts_cl <- counts_split[[cluster]]


# 3. Read pseudo gene score matrix

geneScoreMat <- read.csv(paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/Pseudobulk_GeneScoreMatrix_AllChromosomes.csv"),
					header = TRUE)
rownames(geneScoreMat) <- geneScoreMat[,1]
geneScoreMat[,1] <- NULL
geneScoreMat <- t(geneScoreMat)

# correct rownames
rownames(geneScoreMat) <- gsub("(_ATAC\\.)", "#", rownames(geneScoreMat))
rownames(geneScoreMat) <- gsub("\\.", "-", rownames(geneScoreMat))
tmp <- str_split(rownames(geneScoreMat), pattern = "#", n = 2, simplify=TRUE)
rownames(geneScoreMat) <- paste0(tmp[,2], "__", tmp[,1])

# 3a. Split count data to list with entry per celltype
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(geneScoreMat), 
                                    pattern = "__",  
                                    n = 2), 
                 `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
genescore_split <- split.data.frame(geneScoreMat, 
                                 factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 str_split(rownames(u), '__', simplify = TRUE)[,2]))
# Explore the different components of list
str(genescore_split)

# Subset the genescore matrix to only the specific celltype
#genescore <- genescore_split[[cluster]]


# 4. Get overall correlation across celltypes

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

# correlation for actual data (without randomzation)
rowCor_mat <- sapply(1:ncol(counts), function(i) cor(counts[,i], geneScoreMat[,i], method = "pearson"))
quantile(rowCor_mat, na.rm=TRUE)

# randomize columns in gene score matrix
set.seed(8)
geneScoreMat_rand <- geneScoreMat[sample(nrow(geneScoreMat)),]

# correlation for randomized data
rowCor_mat_rand <- sapply(1:ncol(counts), function(i) cor(counts[,i], geneScoreMat_rand[,i], method = "pearson"))
quantile(rowCor_mat_rand, na.rm=TRUE)

# # randomize multiple times with different seed values
# list_rand <- list()
# for (s in 1:50){
# 	set.seed(s)
# 	# randomize columns in gene score matrix
# 	geneScoreMat_rand <- geneScoreMat[sample(nrow(geneScoreMat)),]

# 	# correlation for randomized data
# 	rowCor_mat_rand <- sapply(1:ncol(counts), function(i) cor(counts[,i], geneScoreMat_rand[,i], method = "pearson"))
# 	quantile(rowCor_mat_rand, na.rm=TRUE)
# 	list_rand[[s]] <- rowCor_mat_rand
# }
# mat_rand <- do.call(rbind, list_rand)
# rowCor_mat_rand_mean <- colMeans(mat_rand)

# 4a. Plot correlation values
plot_df <- data.frame(gene_id = rep(colnames(counts),2),
					corr = c(rowCor_mat, rowCor_mat_rand),
					group = c(rep("actual", ncol(counts)), 
							rep("random", ncol(counts))))
plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
# dataframe for mean value per group
mu <- plot_df %>%
    group_by(group) %>%
    dplyr::summarize(group_mean = mean(corr, na.rm=TRUE))
# perform t-test (normally distributed?) between distributions
tres <- t.test(rowCor_mat, rowCor_mat_rand, 
			alternative="greater")

# all genes, colored according to actual/random data
ggplot(plot_df, aes(x=corr, fill=group, color=group))+
	geom_histogram(position="identity", alpha=0.5) +
	geom_vline(data=mu, aes(xintercept=group_mean, color=group),
             linetype="dashed") +
			 ggtitle(paste(ct, ": all genes, t-test pvalue: ", tres$p.value)) +
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
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_Random_Histogram.pdf"),
	width= 8,
	height=6)

# plot for DE genes only
plot_df <- plot_df[plot_df$de_hit == TRUE,]
mu <- plot_df %>%
    group_by(group) %>%
    dplyr::summarize(group_mean = mean(corr, na.rm=TRUE))
# perform t-test (normally distributed?) between distributions
tres <- t.test(plot_df$corr[plot_df$group == "actual"], 
			plot_df$corr[plot_df$group == "random"], 
			alternative="greater")

ggplot(plot_df, aes(x=corr, fill=group, color=group))+
	geom_histogram(position="identity", alpha=0.5) +
	geom_vline(data=mu, aes(xintercept=group_mean, color=group),
             linetype="dashed") +
			 ggtitle(paste(ct, ": DE hits, t-test pvalue: ", tres$p.value)) +
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
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_Random_HistogramDEhits.pdf"),
	width= 8,
	height=6)

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
ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_MeanVariance.pdf"),
	width= 8,
	height=6)


# 5. Get correlation per celltype
for (ct in celltypes){

	if (ct %in% names(genescore_split)){

	counts_ct <- t(counts_split[[ct]])
	genescore_ct <- t(genescore_split[[ct]])

	# reorder gene score matrix according to ordering of gene expression matrix
	# and subset count matrix for common pseudosamples
	indices <- intersect(rownames(counts_ct), rownames(genescore_ct))
	genescore_ct <- genescore_ct[indices,]
	counts_ct <- counts_ct[indices,]

	# subset both matrices to common set of genes
	colnames(genescore_ct) <- mapIds(org.Hs.eg.db,
						keys=colnames(genescore_ct),
						column="ENSEMBL",
						keytype="SYMBOL",
						multiVals="first")
	indices_col <- intersect(colnames(counts_ct), colnames(genescore_ct))
	genescore_ct <- genescore_ct[,indices_col]
	counts_ct <- counts_ct[,indices_col]

	# correlation for actual data (without randomzation)
	rowCor_mat <- sapply(1:ncol(counts_ct), function(i) cor(counts_ct[,i], genescore_ct[,i], method = "pearson"))

    # randomize columns in gene score matrix
	set.seed(8)
    genescore_ct_rand <- genescore_ct[sample(nrow(genescore_ct)),]
	# correlation for randomized gene score matrix
	rowCor_mat_rand <- sapply(1:ncol(counts_ct), function(i) cor(counts_ct[,i], genescore_ct_rand[,i], method = "pearson"))

	# # randomize multiple times with different seed values
	# list_rand <- list()
	# for (s in 1:50){
	# 	set.seed(s)
	# 	# randomize columns in gene score matrix
	# 	genescore_ct_rand <- genescore_ct[sample(nrow(genescore_ct)),]

	# 	# correlation for randomized data
	# 	rowCor_mat_rand <- sapply(1:ncol(counts_ct), function(i) cor(counts_ct[,i], genescore_ct_rand[,i], method = "pearson"))
	# 	quantile(rowCor_mat_rand, na.rm=TRUE)
	# 	list_rand[[s]] <- rowCor_mat_rand
	# }
	# mat_rand <- do.call(rbind, list_rand)
	# rowCor_mat_rand_mean <- colMeans(mat_rand)


	# Plot correlation values
	plot_df <- data.frame(gene_id = rep(colnames(counts_ct),2),
						corr = c(rowCor_mat, rowCor_mat_rand),
						group = c(rep("actual", ncol(counts_ct)), 
								rep("random", ncol(counts_ct))))
	plot_df$de_hit <- ifelse(plot_df$gene_id %in% df$gene, TRUE, FALSE)
	# dataframe for mean value per group
	mu <- plot_df %>%
		group_by(group) %>%
		dplyr::summarize(group_mean = mean(corr, na.rm=TRUE))
	# perform t-test (normally distributed?) between distributions
	tres <- t.test(rowCor_mat, rowCor_mat_rand, 
				alternative="greater")

	# all genes, colored according to actual/random data
	ggplot(plot_df, aes(x=corr, fill=group, color=group))+
		geom_histogram(position="identity", alpha=0.5) +
		geom_vline(data=mu, aes(xintercept=group_mean, color=group),
				linetype="dashed") +
				ggtitle(paste(ct, ": all genes, t-test pvalue: ", tres$p.value)) +
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
	ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_Random_Histogram_",ct,".pdf"),
		width= 8,
		height=6)

	# plot for DE genes only
	plot_df <- plot_df[plot_df$de_hit == TRUE,]
	mu <- plot_df %>%
		group_by(group) %>%
		dplyr::summarize(group_mean = mean(corr, na.rm=TRUE))
	# perform t-test (normally distributed?) between distributions
	tres <- t.test(plot_df$corr[plot_df$group == "actual"], 
				plot_df$corr[plot_df$group == "random"], 
				alternative="greater")
	ggplot(plot_df, aes(x=corr, fill=group, color=group))+
		geom_histogram(position="identity", alpha=0.5) +
		geom_vline(data=mu, aes(xintercept=group_mean, color=group),
				linetype="dashed") +
				ggtitle(paste(ct, ": DE hits, t-test pvalue: ", tres$p.value)) +
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
	ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_Random_HistogramDEhits_",ct,".pdf"),
		width= 8,
		height=6)

	}

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
	ggsave(filename = paste0(DEA_dir, "plots/03_3b_CorrelationGeneScoreAndGeneExp_Pseudobulk_Random_OneSeed/03_3b_CorrelationGeneScoresAndGeneExp_Pseudobulk_MeanVariance", ct, ".pdf"),
		width= 8,
		height=6)

}



