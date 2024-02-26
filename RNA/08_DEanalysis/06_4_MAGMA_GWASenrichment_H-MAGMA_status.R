##################################################
## Project: SingleNuclei Postmortem Brain
## Date: 17.02.2023
## Author: Nathalie
##################################################
# Compare celltype specific DEGs with previously identified GWAS risk genes 
# of psychiatric disorders (https://pubmed.ncbi.nlm.nih.gov/32152537/)
# with the MAGMA enrichment method

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


## 2. Read in DE genes
# read all genes to have them as background lists and use first n genes
p <- paste0("X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
files <- files[!str_detect(files, "metadata")]
# read in the diff genes of different cell types
de_list <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_list[[ct]] <- de[order(de$padj),]
}
celltypes <- names(de_list)

# lists containing only the gene ids
gene_list <- lapply(de_list, function(x) x$gene[1:250])


## 3. Prepare MAGMA input for each trait and celltype 
# # --> doesn't work, as this file is not accepted by MAGMA
# for (trait in traits){
#   
#   # file path to H-MAGMA results
#   f_trait <- paste0(basedir, "/tables/external/H-MAGMA/",
#                     "Adult_brain_", trait, ".genes.raw")
#   
#   ## Read GWAS/H-MAGMA data
#   no_col <- max(count.fields(f_trait, sep = " ")) # max number of columns per line
#   data <- read.csv(f_trait, skip = 2,
#                    header = FALSE
#                    #col.names=c("gene", "chr", "start", "end", 1:no_col)
#                    ) 
#   data$gene <- word(data$V1, 1)
#   
#   for (ct in celltypes){
#     # background in specific celltype
#     background_ct <- de_list[[ct]]$gene
#     
#     # subset to background set
#     data_ct <- data %>%
#       filter(gene %in% background_ct) %>%
#       dplyr::select(V1)
#     
#     write.table(data_ct, file = paste0(basedir, "/tables/GWASenrichment/MAGMAinput/",
#                                        trait, "_", ct, "_H-MAGMA_filteredTable.genes.raw"),
#                 quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#     
#   }
# }


# 4. Create set-annot file for each cell type
for (ct in celltypes){
  
  # prepare gene set file for this cell type
  gene_set <- data.frame("geneset"=ct,
                         "genes" = paste(gene_list[[ct]], collapse=" "))
  
  write.table(gene_set, file = paste0(basedir, "/tables/GWASenrichment/MAGMAinput/status/",
                                      ct, "_geneset.annot"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # prepare file with continuous variable and all genes
  gene_df <- de_list[[ct]] %>%
    dplyr::select(gene, padj, log2FoldChange) 
  gene_df$padj[is.na(gene_df$padj)] <- 1
  gene_df <- gene_df %>% 
    mutate(log_padj = -log10(padj)) %>%
    mutate(score = log_padj*log2FoldChange)
  write.table(gene_df, file = paste0(basedir, "/tables/GWASenrichment/MAGMAinput/status/",
                                      ct, "_gene_covar.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
}
