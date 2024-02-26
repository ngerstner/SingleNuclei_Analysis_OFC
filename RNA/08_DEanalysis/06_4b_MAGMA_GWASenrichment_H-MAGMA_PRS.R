##################################################
## Project: SingleNuclei Postmortem Brain
## Date: 22.02.2023
## Author: Nathalie
##################################################
# Analyse self-calculated H-MAGMA results

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(fgsea)
library(biomaRt)
library(ggplot2)

## 1. Setup

# Define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"

# GWAS of interest from H-MAGMA results
traits <- c("BIP2021", "MDD", "SCZ2022")

# PRS trait
# prs_trait <- "crossDisorder2019"
# prs_trait <- "MDD"
# prs_trait <- "SCZ2022"
# prs_trait <- "BIP2021"
prs_trait <- "height"

# 1. Read in calculated H-MAGMA results
magma_list <- list()
for (trait in traits){
  
  # find all files for one trait with file pattern
  p <- paste0(trait, "_.*_magmaEnrichment.gsa.out")
  files <- list.files(path = paste0(basedir, "tables/GWASenrichment/MAGMAoutput/", 
                                    prs_trait, "/"),
                      pattern = p,
                      full.names = FALSE)
  # read in the files
  magma_score_list <- list()
  for (f in files){
    m <- read.table(paste0(basedir, "tables/GWASenrichment/MAGMAoutput/",
                           prs_trait, "/", f),
                    skip = 3, header = TRUE)
    ct <- sub(paste0(trait,"_(.*)_magmaEnrichment.gsa.out"),"\\1",f)
    magma_score_list[[ct]] <- m[m$VARIABLE == "score",]
  }
  celltypes <- names(magma_score_list)
  
  magma_score_df <- bind_rows(magma_score_list, .id="celltype")
  # magma_score_df$FDR <- p.adjust(magma_score_df$P, method = "fdr")
  
  magma_list[[trait]] <- magma_score_df
}

magma_df <- bind_rows(magma_list, .id="trait")
magma_df$FDR <- p.adjust(magma_df$P, method="fdr")


## 2. Plot MAGMA p-values as heatmap from all traits and cell types
# magma_df$sig <- ifelse(magma_df$FDR <= 0.1, "*", "")
# magma_df$sig[magma_df$FDR <= 0.05] <- "**"
magma_df$sig[magma_df$FDR <= 0.05] <- "*"
ggplot(magma_df, aes(x = trait, y = celltype, fill = -log10(FDR), label = sig)) +
  geom_tile() +
  geom_text() +
  xlab("GWAS trait") +
  ylab("Celltype") +
  scale_x_discrete(labels = c("Bipolar Disorder", "MDD", "Schizophrenia")) +
  scale_fill_gradient(name = "-log10(FDR)",
                      low="lightgrey", high="#D45E60",
                      limits = c(0,max(1,max(-log10(magma_df$FDR))))) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/06_4_MAGMA_GWASenrichment_H-MAGMA/", 
                         "/06_4_MAGMA_GWASenrichment_H-MAGMA_DE_", 
                         prs_trait, ".pdf"),
       width = 8,
       height = 6)
