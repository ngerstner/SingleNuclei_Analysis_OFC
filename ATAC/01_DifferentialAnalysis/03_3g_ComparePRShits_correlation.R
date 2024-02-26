## Comparison of PRS hits between PRS traits based
## on the correlation between effect sizes

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(tidyr)
library(forcats)
library(viridis)


## 1. Setup

# define pathes to data and folders
basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# reference trait to compare to 
# ref_trait <- "height"
# ref_trait <- "crossDisorder2019"
# ref_trait <- "MDD"
ref_trait <- "SCZ2022"

# PRS trait
# prs_trait <- "crossDisorder2019"
# prs_trait <- "MDD"
# prs_trait <- "SCZ2022"
prs_trait <- "BIP2021"

# outlier removed version of results
outlier <- "_OutlierRemoved"
# outlier <- ""

# gene score matrix normalized?
norm <- "_Normalized"
# norm <- ""

celltypes <- c("Astro_FB", "Astro_PP", "Endothelial", "Exc_L2-3", "Exc_L3-5",
		"Exc_L4-6_1", "Exc_L4-6_2", "In_PVALB_Ba", "In_PVALB_Ch", 
		"In_RELN", "In_SST", "In_VIP", "Microglia", "Oligodendrocyte", "OPC")

## 1. Read in DA genes
# 1.1 Height
res_list <- list()
for (ct in celltypes){
    p <- paste0("GeneScore_WaldTest", outlier, norm, "_", ct, "_", ref_trait, ".csv")
    # read in results of differential gene score matrix analysis
    res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/", p))
    res_list[[ct]] <- res
}
res_list <- res_list[sapply(res_list, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# 2.2 Comparison trait
comp_list <- list()
for (ct in celltypes){
    p <- paste0("GeneScore_WaldTest", outlier, norm, "_", ct, "_", prs_trait, ".csv")
    # read in results of differential gene score matrix analysis
    res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/", p))
    comp_list[[ct]] <- res
}
comp_list <- comp_list[sapply(comp_list, function(x) dim(x)[1]) > 0] # remove empty dataframes from list


# 3. Correlation of effect sizes
corr_effSize <- function(df_cd, df_scz){
  
  # set gene_symbols as rownames
  rownames(df_cd) <- df_cd$gene_name
  rownames(df_scz) <- df_scz$gene_name
  
  # subset dataframes to common list of genes
  common_genes <- intersect(df_cd$gene, df_scz$gene)
  cor_es <- cor.test(df_cd[common_genes,]$log2FC, 
                     df_scz[common_genes,]$log2FC)
  
  return(c("corr_estimate" = cor_es$estimate,
           "pvalue" = cor_es$p.value))
}

df <- data.frame("ct_height" = character(),
                 "ct_comp" = character(),
                 "corr_estimate" = numeric(),
                 "pvalue" = numeric())
for (ct_cd in names(res_list)){
  for(ct_scz in names(comp_list)) {
    
    test_comb <- corr_effSize(as.data.frame(res_list[[ct_cd]]), 
                            as.data.frame(comp_list[[ct_scz]]))
    
    # new_line <- c(ct_cd, ct_scz,
    #               test_comb["corr_estimate.cor"], test_comb["pvalue"])
    # names(new_line) <- c("ct_cd", "ct_scz", "corr_estimate", "pvalue")
    # df <- rbind(df, new_line)
    df <- df %>% add_row("ct_height" = ct_cd, "ct_comp" = ct_scz,
                         "corr_estimate" = test_comb["corr_estimate.cor"], 
                         "pvalue" = test_comb["pvalue"])
  }
}

max_value <- max(abs(df$corr_estimate))
# df$FDR <- p.adjust(df$pvalue, method = "fdr")
# df$label <- ifelse(df$FDR <= 0.05, "*", "")
ggplot(df, aes(x=ct_height, y=ct_comp, fill=corr_estimate)) +
  geom_tile() +
  # geom_text(aes(label = label, size=2.5)) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white", 
                       midpoint = 0, limit = c(-0.75,0.75), space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab(paste0("Cell type ", ref_trait)) +
  ylab(paste0("Cell type ", prs_trait)) +
  theme_light() +
  theme(
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(
  filename = paste0(DEA_dir, "plots/03_3g_ComparePRShits/03_3g_ComparePRShits_", ref_trait, "_correlation_", 
                    prs_trait, ".pdf"),
  width = 12,
  height = 8
)  