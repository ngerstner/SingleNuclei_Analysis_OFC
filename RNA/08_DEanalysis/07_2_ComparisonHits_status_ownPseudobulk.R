## Compare the cell type specific DE hits for disease status
## with the full pseudobulk hits for disease status

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(ggpubr)

## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')
plotdir <- paste0(basedir, "plots/")

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"


## 2. Read in DE genes
de_genes <- list()
all_genes <- list()

# 2.1 Read genes from case control per celltype
p <- paste0("X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
files <- files[!str_detect(files, "metadata")]
# read in the diff genes of different cell types
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  de[is.na(de)] <- 1
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_genes[[ct]] <- de[de$padj <= 0.1,]
  de_genes[[ct]] <- de_genes[[ct]] %>%
    mutate(up_down = ifelse(log2FoldChange > 0, "up", "down"))
  all_genes[[ct]] <- de
}
celltypes <- names(de_genes)

# 2.2 Read genes from case control pseudobulk analysis
file_pseudo <-
  paste0(
    basedir,
    "tables/DESeq2_Wald/PseudobulkComplete/PseudobulkComplete_~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch_Wald_pcFromCorrectedDataAllCellTypes.csv"
  )
de_pseudo <- read.csv(file_pseudo)

# filter genes according to FDR cutoff of 0.1
de_pseudo[is.na(de_pseudo)] <- 1
de_genes_pseudo <- de_pseudo[de_pseudo$padj <= 0.1,]
de_genes_pseudo <- de_genes_pseudo %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "up", "down"))


## 3. Compare cell type hits to pseudobulk hits

# 3.1 Overlap with pseudobulk per cell type
de_genes <- de_genes[lapply(de_genes,nrow)>0]

for (ct in names(de_genes)){
  
  print(ct)
  
  de_genes[[ct]] <- de_genes[[ct]] %>%
    mutate(de_pseudobulk = ifelse(gene %in% de_genes_pseudo$gene, TRUE, FALSE)) %>%
    left_join(de_pseudo[c("gene", "log2FoldChange", "padj")], by="gene") %>%
    mutate(up_down_pseudo = ifelse(log2FoldChange.y > 0, "up", "down")) %>%
    mutate(dir_same = (up_down == up_down_pseudo))
}

df_all <- bind_rows(de_genes, .id = "celltype")
df_all$celltype <- factor(df_all$celltype, levels = celltypes)

# check genes with contradictory directions of regulation in pseudubulk and cell type
df_all[!df_all$dir_same,] # -> none of them is significant in pseudobulk

# count number of genes per group
df_all_plot <- df_all %>%
  group_by(celltype, up_down, de_pseudobulk) %>%
  count()
df_all_plot$up_down <- factor(df_all_plot$up_down, levels = c("up", "down"))

# 3.2 Plot Overlap of pseudobulk and cell type hits
ggplot(df_all_plot, aes(x=celltype, y=n, fill=up_down)) +
  geom_bar(stat = "identity", aes(alpha=de_pseudobulk)) +
  facet_wrap(~up_down, nrow=2,
             labeller = labeller(up_down = c("down"= "Downregulated genes",
                                          "up" = "Upregulated genes"))) +
  scale_alpha_manual("Pseudobulk data", 
                     breaks = c(TRUE, FALSE),
                     values = c(1, 0.5),
                     labels = c("diff. exp.", "not diff. exp.")) +
  scale_fill_manual("Directionality of regulation",
                    breaks = c("up", "down"),
                    values = c("red", "darkblue"),
                    labels = c("upregulated", "downregulated")) +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab("Celltype") +
  ylab("Number of DE genes") +
  theme_light() +
  theme(
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
ggsave(
  filename = paste0(plotdir, "/07_2_ComparisonHits_status_ownPseudobulk/",
                    "/07_2_Barplot_DE_status_ownPseudobulk_", test, ".pdf"),
  width = 12,
  height = 8
)  


## 4. Compare pseudobulk hits to cell type hits
df_all_genes <- bind_rows(all_genes, .id="celltype")
de_genes_pseudo <- de_genes_pseudo %>%
  mutate(de_celltype = ifelse(gene %in% df_all$gene, TRUE, FALSE)) 

de_genes_pseudo_plot <- de_genes_pseudo %>%
  group_by(de_celltype, up_down) %>%
  count() %>%
  mutate(data = "pseudobulk")
de_genes_pseudo_plot$up_down <- factor(de_genes_pseudo_plot$up_down, levels = c("up", "down"))

# 3.2 Plot Overlap of pseudobulk and cell type hits
ggplot(de_genes_pseudo_plot, aes(x = data, y=n, fill=up_down)) +
  geom_bar(stat = "identity", aes(alpha=de_celltype)) +
  facet_wrap(~up_down, nrow=2,
             labeller = labeller(up_down = c("down"= "Downregulated genes",
                                             "up" = "Upregulated genes"))) +
  scale_alpha_manual("Celltype specific hit", 
                     breaks = c(TRUE, FALSE),
                     values = c(1, 0.5),
                     labels = c("diff. exp.", "not diff. exp.")) +
  scale_fill_manual("Directionality of regulation in pseudobulk",
                    breaks = c("up", "down"),
                    values = c("red", "darkblue"),
                    labels = c("upregulated", "downregulated")) +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  ylab("Number of DE genes") +
  theme_light() +
  theme(
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
ggsave(
  filename = paste0(plotdir, "/07_2_ComparisonHits_status_ownPseudobulk/",
                    "/07_2_Barplot_DE_ownPseudobulk_status_", test, ".pdf"),
  width = 10,
  height = 8
)  


