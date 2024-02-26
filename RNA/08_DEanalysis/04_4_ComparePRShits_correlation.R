## Comparison of PRS results between traits looking
## at the correlation between effect sizes

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(Hmisc)
library(corrplot)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

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

## 1. Read in DE genes
# 1.1 Height
p <- paste0("PMI_Wald_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
            ref_trait, ".csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_Wald/PRSext"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_height <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_Wald/PRSext/",f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_height[[ct]] <- de
}
de_height <- de_height[sapply(de_height, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# 2.2 Comparison trait
p <- paste0("PMI_Wald_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
            prs_trait, ".csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_Wald/PRSext"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_comp <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_Wald/PRSext/",f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_comp[[ct]] <- de
}
celltypes <- names(de_comp)
de_comp <- de_comp[sapply(de_comp, function(x) dim(x)[1]) > 0] # remove empty dataframes from list


# 3. Correlation of effect sizes
corr_effSize <- function(df_cd, df_scz){
  
  # set gene_symbols as rownames
  rownames(df_cd) <- df_cd$gene
  rownames(df_scz) <- df_scz$gene
  
  # subset dataframes to common list of genes
  common_genes <- intersect(df_cd$gene, df_scz$gene)
  cor_es <- cor.test(df_cd[common_genes,]$log2FoldChange, 
                     df_scz[common_genes,]$log2FoldChange)
  
  return(c("corr_estimate" = cor_es$estimate,
           "pvalue" = cor_es$p.value))
}

df <- data.frame("ct_height" = character(),
                 "ct_comp" = character(),
                 "corr_estimate" = numeric(),
                 "pvalue" = numeric())
for (ct_cd in names(de_height)){
  for(ct_scz in names(de_comp)) {
    
    test_comb <- corr_effSize(de_height[[ct_cd]], de_comp[[ct_scz]])
    
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
  filename = paste0(basedir, "plots/04_4_ComparePRShits/04_4_ComparePRShits_", ref_trait, "_correlation_", 
                    prs_trait, ".pdf"),
  width = 12,
  height = 8
)  