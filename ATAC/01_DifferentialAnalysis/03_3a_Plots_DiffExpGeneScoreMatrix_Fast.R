## Plots for DA results on gene score level between cases and controls

# run with conda environment sc-atac

library(dplyr)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# outlier removed version of results
outlier <- "_OutlierRemoved"
# outlier <- ""

# gene score matrix normalized?
norm <- "_Normalized"
# norm <- ""

celltypes <- c("Astro_FB", "Astro_PP", "Endothelial", "Exc_L2-3", "Exc_L3-5",
		"Exc_L4-6_1", "Exc_L4-6_2", "In_RELN", "In_PVALB_Ba", "In_PVALB_Ch", 
		"In_SST", "In_VIP", "Microglia", "Oligodendrocyte", "OPC")

# 1. Read gene score results
res_list <- list()
for (ct in celltypes){
	p <- paste0("GeneScore_WaldTest", outlier, norm, "_", ct, ".csv")
	# read in results of differential gene score matrix analysis
	res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/", p))
	res_list[[ct]] <- res
}
res_df <- bind_rows(res_list, .id="celltype")

# 2. Read DE results
p <- "X6.Batch_Wald_pcFromCorrectedDataAllCellTypes.csv"
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
#de_ct <- de_ct[sapply(de_ct, function(x) dim(x)[1]) > 0] # remove celltypes with zero hits
# list to dataframe
df <- bind_rows(de_ct, .id="celltype")
df$celltype <- factor(df$celltype, levels = celltypes_de)
# subset to the significant hits
df$padj[is.na(df$padj)] <- 1
df_sig <- df[df$padj <= 0.1,] # TODO: --> FIX NAs here
head(df_sig)

# 3. Add FDR values per celltype

fdr_perCelltype <- function(data_df, pvalue_column){

	data_list <- split(data_df, f = data_df$celltype)
	# add column with corrected p-values for each celltype
	data_list <- lapply(data_list, function(x) x %>% dplyr::mutate(FDR = p.adjust(x[[pvalue_column]], method="fdr")))
	out_df <- bind_rows(data_list)
	out_df$celltype <- factor(out_df$celltype, levels=celltypes)
	out_df$sig <- ifelse(out_df$FDR <= 0.1, TRUE, FALSE)

	return(out_df)
}

# on merged object with DE hits
df_merge_de <- df_sig %>%
	inner_join(res_df, by=c("celltype", "X"="gene_name"))
de_df <- fdr_perCelltype(df_merge_de, "pvalue.y") 	

de_df_write <- de_df %>%
	filter(sig) %>%
	dplyr::select(celltype, gene, log2FC, statistic, pvalue.y, FDR, X) %>%
	rename(log2FoldChange = log2FC, stat = statistic, pvalue = pvalue.y, padj = FDR, gene_name = X)
write.csv(de_df_write, paste0(DEA_dir, "tables/GeneScoreMatrix/DA-DEhits_allCellTypes_status.csv"))

# for all genes
df_merge_all <- df %>%
	inner_join(res_df, by=c("celltype", "X"="gene_name"))
all_df <- fdr_perCelltype(df_merge_all, "pvalue.y")

all_df_write <- all_df %>%
	filter(sig) %>%
	dplyr::select(celltype, gene, log2FC, statistic, pvalue.y, FDR, X) %>%
	rename(log2FoldChange = log2FC, stat = statistic, pvalue = pvalue.y, padj = FDR, gene_name = X)
write.csv(all_df_write, paste0(DEA_dir, "tables/GeneScoreMatrix/DAhits_allCellTypes_status.csv"))

# for all genes, without subsetting to the ones tested in RNA
df_merge_all <- df %>%
	right_join(res_df, by=c("celltype", "X"="gene_name"))
all_df_unfiltered <- fdr_perCelltype(df_merge_all, "pvalue.y")

all_df_write <- all_df_unfiltered %>%
	filter(sig) %>%
	dplyr::select(celltype, gene, log2FC, statistic, pvalue.y, FDR, X) %>%
	rename(log2FoldChange = log2FC, stat = statistic, pvalue = pvalue.y, padj = FDR, gene_name = X)
write.csv(all_df_write, paste0(DEA_dir, "tables/GeneScoreMatrix/DAhits_allCellTypes_status_unfiltered.csv"))


# 4. Plot results for gene score matrix

plot_genescore_results <- function(data_df, plot_prefix, include_nonSig = TRUE){

	# dot/volcano plot of genes
	if(include_nonSig){
		data_df <- arrange(data_df, sig)
		data_df$orderrank <- rank(data_df$sig,ties.method="first")
		dot_g <- ggplot(data_df, aes(x=celltype, y = log2FC,
			size = -log10(FDR), col = sig)) +
		geom_point() + 
		scale_x_discrete(drop=FALSE,
			guide = guide_axis(angle=45)) +
		scale_size_continuous(name = "-log10(FDR)") +
		scale_color_manual(name = "Diff. accessible",
					breaks = c(TRUE, FALSE),
					labels = c("sig.", "not sig."),
					values = c("orange", "darkgrey")) +
		xlab("Celltype") +
		ylab("log2(Fold Change)") +
		theme_light() +
		theme(
			axis.title.x = element_blank(),
			axis.title.y = element_text(size=13),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12),
			legend.position = "bottom"
		)
		ggsave(filename = paste0(plot_prefix, "_DotPlot_showAllGenes.pdf"),
		plot = dot_g,
		width=8,
		height=6)
	}else{
		data_df_filt <- filter(data_df, sig)
		dot_g <- ggplot(data_df_filt, aes(x=celltype, y = log2FC,
			size = -log10(FDR))) +
		geom_point(col = "orange") + 
		scale_x_discrete(drop=FALSE,
			guide = guide_axis(angle=45)) +
		scale_size_continuous(name = "-log10(FDR)") +
		xlab("Celltype") +
		ylab("log2(Fold Change)") +
		theme_light() +
		theme(
			axis.title.x = element_text(size=13),
			axis.title.y = element_text(size=13),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12),
			legend.position = "bottom"
		)
		ggsave(filename = paste0(plot_prefix, "_DotPlot_showSigGenes.pdf"),
		plot = dot_g,
		width=8,
		height=6)
	}


	# barplot of number of hits per celltype
	n_hits <- data_df %>%
			dplyr::group_by(celltype, sig) %>%
			dplyr::count() %>%
			filter(sig)
	n_hits$celltype <- factor(n_hits$celltype, levels=levels(data_df$celltype))

	bar_g <- ggplot(n_hits, aes(x=celltype, y = n)) +
		geom_bar(stat = "identity", fill = "darkgrey") + 
		scale_x_discrete(drop=FALSE,
			guide = guide_axis(angle=45)) +
		# xlab("Celltype") +
		ylab("Number of sig. genes") +
		theme_light() +
		theme(
			#axis.title.x = element_text(size=15),
			axis.title.y = element_text(size=13),
			#axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			legend.text = element_text(size=12)
		)
	ggsave(filename = paste0(plot_prefix, "_BarplotNumberGenes.pdf"),
		width=8,
		height=6)

	# dot and barplot combined
	plot_grid(bar_g, dot_g, 
			ncol = 1, align = "v",
			rel_heights = c(0.25,0.75))
	if(include_nonSig){
		ggsave2(filename = paste0(plot_prefix, "_DotPlotAndBarplotNumberGenes_showAllGenes.pdf"),
		#plot = comb_g,
		width=8,
		height=8)
	}else{
	ggsave2(filename = paste0(plot_prefix, "_DotPlotAndBarplotNumberGenes_showSigGenes.pdf"),
		#plot = comb_g,
		width=8,
		height=8)
	}

	# dot plot plotting DE p-value against gene score p-value
	ggplot(data_df, aes(x=-log10(pvalue.x), y =-log10(pvalue.y)))+
		geom_point(aes(col=sig)) +
		geom_smooth(method="lm")+
		xlab("P-value (-log10) differential expression") +
		ylab("P-value (-log10) differential gene score") +
		scale_color_manual(name = "Diff. accessible",
					breaks = c(TRUE, FALSE),
					labels = c("sig.", "not sig."),
					values = c("orange", "darkgrey")) +
		theme_light() +
		theme(
			axis.title.x = element_text(size=15),
			axis.title.y = element_text(size=15),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12),
			legend.position = "bottom"
		)
	ggsave(filename = paste0(plot_prefix, "_PvalueComparison.pdf"),
		width= 6,
		height=6)

	# dot plot plotting DE fold changes against gene score fold changes
	ggplot(data_df, aes(x=log2FoldChange, y =log2FC))+
		geom_point(aes(col=sig)) +
		geom_smooth(method="lm")+
		xlab("log2FC differential expression") +
		ylab("log2FC differential gene score") +
		scale_color_manual(name = "Diff. accessible",
					breaks = c(TRUE, FALSE),
					labels = c("sig.", "not sig."),
					values = c("orange", "darkgrey")) +
		theme_light() +
		theme(
			axis.title.x = element_text(size=15),
			axis.title.y = element_text(size=15),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12),
			legend.position = "bottom"
		)
	ggsave(filename = paste0(plot_prefix, "_FoldChangeComparison.pdf"),
		width= 6,
		height=6)



	# histogram of p-values
	ggplot(data_df, aes(x=-log10(pvalue.y)))+
		geom_histogram() +
		xlab("P-values (-log10) differential gene score") +
		theme_light() +
		theme(
			axis.title.x = element_text(size=15),
			axis.title.y = element_text(size=15),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12)
		)
	ggsave(filename = paste0(plot_prefix, "_HistogramGeneScorePvalues.pdf"),
		width= 6,
		height=6)


	ggplot(data_df, aes(x=-log10(pvalue.x)))+
		geom_histogram() +
		xlab("P-values (-log10) differential expression") +
		theme_light() +
		theme(
			axis.title.x = element_text(size=15),
			axis.title.y = element_text(size=15),
			axis.text.x = element_text(size=12),
			axis.text.y = element_text(size=12),
			legend.text = element_text(size=12)
		)
	ggsave(filename = paste0(plot_prefix, "_HistogramExpressionPvalues.pdf"),
		width= 6,
		height=6)

}

# run plot function for gene score referring to a DE hit
plot_genescore_results(de_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_DEhits", outlier, norm), include_nonSig=TRUE)
plot_genescore_results(de_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_DEhits", outlier, norm), include_nonSig=FALSE)


# run plot function for all gene scores
plot_genescore_results(all_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_AllGenes", outlier, norm), include_nonSig=TRUE)
plot_genescore_results(all_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_AllGenes", outlier, norm), include_nonSig=FALSE)

# run plot function for all gene scores
plot_genescore_results(all_df_unfiltered, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_AllGenesUnfiltered", outlier, norm), include_nonSig=TRUE)
plot_genescore_results(all_df_unfiltered, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_AllGenesUnfiltered", outlier, norm), include_nonSig=FALSE)

