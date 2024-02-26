## Plots for DA results on gene score level between PRS risk groups

# run with conda environment sc-atac

library(dplyr)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(org.Hs.eg.db)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# outlier removed version of results
outlier <- "_OutlierRemoved"
# outlier <- ""

# gene score matrix normalized?
norm <- "_Normalized"
# norm <- ""

traits <- c("crossDisorder2019", "BIP2021", "height", "MDD", "SCZ2022")

celltypes <- c("Astro_FB", "Astro_PP", "Endothelial", "Exc_L2-3", "Exc_L3-5",
		"Exc_L4-6_1", "Exc_L4-6_2", "In_PVALB_Ba", "In_PVALB_Ch", 
		"In_RELN", "In_SST", "In_VIP", "Microglia", "Oligodendrocyte", "OPC")

# 1. Read gene score results
res_list <- list()
for (trait in traits){
    res_list[[trait]] <- list()
    for (ct in celltypes){
        p <- paste0("GeneScore_WaldTest", outlier, norm, "_", ct, "_", trait, ".csv")
        # read in results of differential gene score matrix analysis
        res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/", p))
        res_list[[trait]][[ct]] <- res
    }
	res_list[[trait]] <- res_list[[trait]][sapply(res_list[[trait]], 
                                                  function(x) dim(x)[1]) > 0] # remove empty dataframes from list
}
#res_df <- bind_rows(res_list, .id="celltype")

# 1a. List to dataframe
df_onelevel <- lapply(res_list, function(x) bind_rows(x, .id="celltype"))
df <- bind_rows(df_onelevel, .id = "prs_trait")
df$celltype <- factor(df$celltype, levels = celltypes)
df$prs_trait <- factor(df$prs_trait, levels = traits)

# 2. Read DE results
de_ct <- list()
for (prs_trait in traits){
  p <- paste0("sig0.1_Wald_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
              prs_trait, ".csv")
  files <- list.files(path = paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/PRSext"),
                      pattern = p,
                      full.names = FALSE)
  # read in the diff genes of different cell types
  de_ct[[prs_trait]] <- list()
  for (f in files){
    de <- read.csv(paste0(basedir, "scripts/RNA/08_DEanalysis/DESeq2_RELN/tables/DESeq2_Wald/PRSext/",f))
    ct <- sub("(\\w*)_~.*","\\1",f)
    de_ct[[prs_trait]][[ct]] <- de
  }
  de_ct[[prs_trait]] <- de_ct[[prs_trait]][sapply(de_ct[[prs_trait]], 
                                                  function(x) dim(x)[1]) > 0] # remove empty dataframes from list
}

## 3. List to dataframe
de_onelevel <- lapply(de_ct, function(x) bind_rows(x, .id="celltype"))
de <- bind_rows(de_onelevel, .id = "prs_trait")
de$celltype <- factor(de$celltype, levels = celltypes)
de$prs_trait <- factor(de$prs_trait, levels = traits)


# 3. Add FDR values per celltype

fdr_perCelltype <- function(data_df, pvalue_column){

	out_list <- list()
	for (trait in traits){
		# subset df to rows of respective PRS trait
		data_df_trait <- data_df %>% filter(prs_trait == trait)
		
		# split df into a list with df per celltype
		data_list <- split(data_df_trait, f = data_df_trait$celltype)
		data_list <- data_list[sapply(data_list, 
                                                  function(x) dim(x)[1]) > 0] # remove empty dataframes from list
		# add column with corrected p-values for each celltype
		data_list <- lapply(data_list, function(x) x %>% dplyr::mutate(FDR = p.adjust(x[[pvalue_column]], method="fdr")))
		
		# TODO: adapt this
		out_df <- bind_rows(data_list)
		out_df$celltype <- factor(out_df$celltype, levels=celltypes)
		out_df$sig <- ifelse(out_df$FDR <= 0.1, TRUE, FALSE)
		out_list[[trait]] <- out_df
	}
	out_df_alltraits <- bind_rows(out_list)

	return(out_df_alltraits)
}

# on merged object with DE hits
df_merge_de <- df %>%
	inner_join(de, by=c("prs_trait", "celltype", "gene_name"="X"))
de_df <- fdr_perCelltype(df_merge_de, "pvalue.x") 	

de_df_write <- de_df %>%
	filter(sig) %>%
	dplyr::select(prs_trait, celltype, gene, log2FC, statistic, pvalue.x, FDR, gene_name) %>%
	rename(log2FoldChange = log2FC, stat = statistic, pvalue = pvalue.x, padj = FDR)
write.csv(de_df_write, paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/DA-DEhits_allCellTypes_PRSext.csv"))

# for all genes
all_df <- fdr_perCelltype(df, "pvalue")
all_df <- all_df %>%
			left_join(de, by=c("prs_trait", "celltype", "gene_name"="X")) %>%
			mutate(DEhit = ifelse(!is.na(gene), TRUE, FALSE))
all_df$gene <- mapIds(org.Hs.eg.db, keys = all_df$gene_name,
                          keytype="SYMBOL", column = "ENSEMBL")

all_df_write <- all_df %>%
	filter(sig) %>%
	dplyr::select(prs_trait, celltype, gene, log2FC, statistic, pvalue.x, FDR, gene_name) %>%
	dplyr::rename(log2FoldChange = log2FC, stat = statistic, pvalue = pvalue.x, padj = FDR)
write.csv(all_df_write, paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/DAhits_allCellTypes_PRSext.csv"))


# 4. Plot results for gene score matrix

plot_genescore_results_sigHits <- function(data_df, plot_prefix){

	data_df <- data_df %>% 
				mutate(prs_trait = factor(prs_trait, 
									levels = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021")))
	prs_names <- c(
                    `crossDisorder2019` = "Cross-Disorder",
                    `SCZ2022` = "SCZ",
					`height` = "Height",
                    `MDD` = "MDD",
                    `BIP2021` = "BIP"
                    )

	# dot/volcano plot with facet wrap
	# color corresponds to PRS trait
	data_df_sig <- data_df %>%
					filter(sig == TRUE)
	ggplot(data_df_sig, aes(x=celltype, y=log2FC, 
               size = -log10(FDR), color=prs_trait)) +
		geom_point() +
		scale_x_discrete(drop=FALSE,
						guide = guide_axis(angle = 90)) +
		scale_size_continuous(name = "-log10(FDR)") +
		scale_color_manual(name = "GWAS Study",
                     breaks = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021"),
                     labels = c("Cross-Disorder", "Schizophrenia", "Height", "Major Depressive Disorder", "Bipolar Disorder"),
                     values = c("#473335", "#548687", "#CFCFA2", "#B0413E", "#FCAA67"),
                    guide = "none") +
		facet_wrap(~prs_trait, drop = FALSE) +
		xlab("Cell type") +
		ylab("log2(Fold Change)") +
		theme_light() +
		theme(
			strip.text.x = element_text(size = 15),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size = 15),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12)
		)
	ggsave(
	filename = paste0(
		plot_prefix, "_DotPlot_sigHits.pdf"
	),
	width = 12,
	height = 8
	)

	df_nr_all <- data_df_sig %>%
		dplyr::count(celltype, prs_trait)
	ggplot(df_nr_all, aes(x=celltype, y=n, fill=prs_trait)) +
		geom_col(position=position_dodge2(preserve = "single")) +
		scale_y_log10() +
		scale_x_discrete(drop=FALSE,
						guide = guide_axis(angle = 90)) +
		scale_fill_manual(name = "GWAS Study",
                     breaks = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021"),
                     labels = c("Cross-Disorder", "Schizophrenia", "Height", "Major Depressive Disorder", "Bipolar Disorder"),
                     values = c("#473335", "#548687", "#CFCFA2", "#B0413E", "#FCAA67"),
                    guide = "none") +
		facet_wrap(~prs_trait, labeller = as_labeller(prs_names)) +#, scales = "free_y") +
		xlab("Cell type") +
		ylab("Number of DA genes (FDR 0.1)") +
		theme_light() +
		theme(
			strip.text.x = element_text(size = 15),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size = 15),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12)
		)
	ggsave(
		filename = paste0(
		plot_prefix, "_Barplot-NumberDAgenes_sigHits.pdf"
		),
		width = 12,
		height = 8
	)

}

plot_genescore_results_sigHits_DEcolor <- function(data_df, plot_prefix){

	data_df <- data_df %>% 
				mutate(prs_trait = factor(prs_trait, 
									levels = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021")))
	prs_names <- c(
                    `crossDisorder2019` = "Cross-Disorder",
                    `SCZ2022` = "Schizophrenia",
                    `MDD` = "Major Depressive Disorder",
                    `BIP2021` = "Bipolar Disorder",
					`height` = "Height"
                    )

	# dot/volcano plot with facet wrap
	data_df_sig <- data_df %>%
					filter(sig == TRUE)
	data_df_sig_DE <- data_df_sig %>%
					filter(DEhit == TRUE) 
	ggplot(data_df_sig, aes(x=celltype, y=log2FC, 
               size = -log10(FDR))) +
		geom_point(aes(color=prs_trait)) +
		geom_point(data = data_df_sig_DE, mapping = aes(x=celltype, y=log2FC,
						size = -log10(FDR)), color="grey") +
		scale_x_discrete(drop=FALSE,
						guide = guide_axis(angle = 90)) +
		scale_size_continuous(name = "-log10(FDR)") +
		scale_color_manual(name = "GWAS Study",
                     breaks = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021"),
                     labels = c("Cross-Disorder", "Schizophrenia", "Height", "Major Depressive Disorder", "Bipolar Disorder"),
                     values = c("#473335", "#548687", "#CFCFA2", "#B0413E", "#FCAA67"),
                    guide = "none") +
		facet_wrap(~prs_trait, labeller = as_labeller(prs_names)) +
		xlab("Cell type") +
		ylab("log2(Fold Change)") +
		theme_light() +
		theme(
			strip.text.x = element_text(size = 15),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size = 15),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12)
		)
	ggsave(
	filename = paste0(
		plot_prefix, "_DotPlot_sigHits_DEcolor.pdf"
	),
	width = 12,
	height = 8
	)

}


# run plot function for gene score referring to a DE hit
plot_genescore_results_sigHits(de_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_PRSext_DEhits", outlier, norm))

# run plot function for all gene scores
plot_genescore_results_sigHits(all_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_PRSext_AllGenes", outlier, norm))

# plot for all genes with color indicating gene is a DE hit
plot_genescore_results_sigHits_DEcolor(all_df, paste0(DEA_dir, "plots/03_3a_DifferentialGeneScores_PRSext_AllGenes", outlier, norm))


