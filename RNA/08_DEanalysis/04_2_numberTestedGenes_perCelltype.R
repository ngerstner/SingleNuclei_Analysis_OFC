## Plot the number of tested genes per cell type
## as part of the downsampling analysis

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(viridis)
library(pheatmap)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')


## 2. Read in DE genes

# 2.1 Read subsampled files for different percentiles
de_ct <- list()
for (perc in c("perc25","perc50","perc75")){
  # list files with DE results in celltypes
  p <- paste0("X6.Batch_Wald_pcFromCorrectedDataAllCellTypes_",perc,".csv")
  files <- list.files(path = file.path(basedir, "tables/DESeq2_Wald/downsampled/"),
                      pattern = p,
                      full.names = FALSE)
  # read in the diff genes of different cell types
  de_ct[[perc]] <- list()
  for (f in files){
    de <- read.csv(file.path(basedir, "tables/DESeq2_Wald/downsampled",f))
    ct <- sub("(\\w*)_~.*","\\1",f)
    de_ct[[perc]][[ct]] <- de[c("gene", "padj")]
  }
}


# 2.2 Read files with all cells - not subsampled
# list files with DE results in celltypes
p <- "X6.Batch_Wald_pcFromCorrectedDataAllCellTypes.csv"
files <- list.files(path = file.path(basedir, "tables/DESeq2_Wald/"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_ct[["perc100"]] <- list()
genes_ct <- list()
for (f in files){
  de <- read.csv(file.path(basedir, "tables/DESeq2_Wald",f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_ct[["perc100"]][[ct]] <- de[c("gene", "padj")]
  genes_ct[[ct]] <- de$gene
}
celltypes <- names(de_ct[["perc100"]])

## 3. Get number of tested genes per celltype and percentile
n_ct_onelevel <- lapply(de_ct, function(x) stack(sapply(x, nrow)))
n_ct <- bind_rows(n_ct_onelevel,.id = "percentile")
n_ct$percentile <- factor(n_ct$percentile, 
                                 levels = c("perc25", "perc50", "perc75", "perc100"))


## 4. Plot number of tested genes for different percentiles
ggplot(n_ct, aes(x=ind, y=values, fill=percentile)) +
  # geom_bar(position="dodge", stat="identity") +
  geom_col(position=position_dodge2(preserve = "single")) +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  scale_fill_viridis(discrete = T,
                     name="Downsampling",
                     breaks = c("perc25", "perc50", "perc75", "perc100"),
                     labels = c("25%", "50%", "75%", "100%")) +
  xlab("Cell type") +
  ylab("Number of tested genes in DE analysis") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_2_NumberOfTestedGenes.pdf"),
       width = 8,
       height = 6)


## 5. Plot overlap of genes between celltypes using all cells
overlap_genes <- crossprod(table(stack(genes_ct)))

pdf(file = paste0(basedir, "plots/04_2_Overlap_TestedGenes.pdf"),
    width = 10, height = 8)
pheatmap(overlap_genes,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

pdf(file = paste0(basedir, "plots/04_2_Overlap_TestedGenes_clust.pdf"),
    width = 10, height = 8)
pheatmap(overlap_genes)
dev.off()

## 6. Create plot of tested genes only for 100 percent
tested_genes <- de_ct[["perc100"]]
number_tested <- stack(sapply(tested_genes, nrow))
ggplot(number_tested, aes(x=ind, y=values)) +
  geom_bar(stat="identity", fill="steelblue") +
  #geom_col(position=position_dodge2(preserve = "single")) +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab("Cell type") +
  ylab("Number of tested genes in DE analysis") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_2_NumberOfTestedGenes_DEanalysis.pdf"),
       width = 8,
       height = 6)
