## Plot the FDR + fold change for some genes across cell types
## to see if dysregulation is cell type specific 

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(pheatmap)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"


## 2. Read in DE genes
# read all genes to have them as background lists and use first n genes
p <- paste0("X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_up <- list()
de_down <- list()
genes <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de$padj[is.na(de$padj)] <- 1
  genes[[ct]] <- de
  de_up[[ct]] <- de[de$padj <= 0.1 & de$log2FoldChange > 0,]
  de_down[[ct]] <- de[de$padj <= 0.1 & de$log2FoldChange < 0,]
}
celltypes <- names(de_up)


# 3. Find gene up and downregulated in different cell types
df_up <- bind_rows(de_up, .id="celltype")
df_down <- bind_rows(de_down, .id="celltype")

up_down_genes <- intersect(df_up$gene, df_down$gene)
# -> no gene that is sig. up and downregulated in different cell types


# 4. Find gene with cell-type specific dysregulation
df_hits <- bind_rows(list("up"=df_up, "down"=df_down), .id="regulation")
df_genes <- bind_rows(genes, .id="celltype") %>%
  filter(gene %in% df_hits$gene) %>%
  dplyr::select(gene, log2FoldChange, padj, celltype) 

df_hits_ext <- df_hits %>%
  left_join(df_genes, by="gene")



## Plots of SLIT2 expression level
ens_slit2 <- "ENSG00000145147"

# Read count and metadata
counts <- readRDS(paste0(basedir, "data/counts.rds"))
metadata <- readRDS(paste0(basedir, "data/metadata.rds"))
metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))

# Split count data to list with entry per celltype
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

slit2 <- lapply(counts_split, function(x) x[ens_slit2,])
slit2_df <- bind_rows(slit2, .id="celltype")
slit2_df_long <- pivot_longer(slit2_df,
                              cols = starts_with("SU"),
                              names_to = "Sample.ID",
                              values_to = "Reads")
slit2_df_long <- left_join(slit2_df_long, metadata[c("sample_id", "Sex", "Status")] %>% distinct(),
                           by=c("Sample.ID"="sample_id"))

# plot expression levels across all cell types
ggplot(slit2_df_long, aes(x=celltype, y=Reads, color=Sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=Status),
             position = position_jitterdodge())

# plot expression level only in Exc_L4-6_1
slit2_df_ct <- slit2_df_long %>%
  filter(celltype == "Exc_L4-6_1")
ggplot(slit2_df_ct, aes(x=Sex, y=Reads, fill=Status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge())


# Plot DE results of gene across all cell types

plot_fc_genes <- function(gene_id, gene_name){
  ens_gene <- gene_id
  
  ## Plot fold change of gene across all cell types
  d_list <- list()
  for (c in celltypes){
    d <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/",
                         c, "_~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch_", 
                         test, "_pcFromCorrectedDataAllCellTypes.csv")
    ) %>%
      dplyr::filter(gene == ens_gene)
    d_list[[c]] <- d
  }
  d_df <- bind_rows(d_list, .id = "celltype")
  d_df$celltype <- factor(d_df$celltype, levels = celltypes)
  
  ggplot(d_df, aes(x=celltype, y=log2FoldChange, size=-log10(padj))) +
    geom_point() +
    scale_x_discrete(name = "Cell type",
                     drop=FALSE,
                     guide = guide_axis(angle = 45)) +
    scale_size_continuous(name = "-log10(FDR)") +
    theme_light() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      strip.text.x = element_text(size = 12)
    )
  ggsave(filename = paste0(basedir, "plots/09_", gene_name, "_FC.pdf"),
         width = 8,
         height = 6)
}

plot_fc_genes("ENSG00000145147", "SLIT2")
plot_fc_genes("ENSG00000152214", "RIT2")
plot_fc_genes("ENSG00000152583", "SPARCL1")
plot_fc_genes("ENSG00000035862", "TIMP2")
plot_fc_genes("ENSG00000028277", "POU2F2")
plot_fc_genes("ENSG00000141741", "MIEN1")
plot_fc_genes("ENSG00000185046", "ANKS1B")
plot_fc_genes("ENSG00000096060", "FKBP5")
plot_fc_genes("ENSG00000102974", "CTCF")
plot_fc_genes("ENSG00000102547", "CAB39L")
plot_fc_genes("ENSG00000225746", "MEG8")
plot_fc_genes("ENSG00000184156", "KCNQ3")
