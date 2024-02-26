## Plot DE hits between cases and controls

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"


## 2. Read in DE genes
# find genes
# p <- "noiseLeft.csv"
# p <- "pcFromCorrectedData.csv"
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_ct <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_ct[[ct]] <- de
}
celltypes <- names(de_ct)
de_ct <- de_ct[sapply(de_ct, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# list to dataframe
df <- bind_rows(de_ct, .id = "celltype")
df$celltype <- factor(df$celltype, levels = celltypes)
write.csv(df, paste0(basedir, "tables/DESeq2_", test, "/",
                    "04_DEhits_allCellTypes_status.csv"))

## 3. Plot DE genes

## plot of all DE genes
ggplot(df, aes(x=celltype, y=log2FoldChange, 
               size = -log10(padj))) +
  geom_point() +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  scale_size_continuous(name = "-log10(FDR)") +
  xlab("Cell type") +
  ylab("log2(Fold Change)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_1_DEgenes_FDR0.1_", test, "_pcFromCorrectedDataAllCellTypes.pdf"),
       width = 8,
       height = 6)

## barplot number of DE genes
df_nr_all <- df %>%
  dplyr::count(celltype)
df_nr_all$celltype <- factor(df_nr_all$celltype, levels = celltypes)
ggplot(df_nr_all, aes(x=celltype, y=n)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab("Cell type") +
  ylab("Number of DE genes (FDR 0.1)") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_1_NumberDEgenes_FDR0.1_", test, "_pcFromCorrectedDataAllCellTypes.pdf"),
      width = 10,
      height = 8)

## barplot number of DE genes separated up and downregulation
df_nr <- df %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "upregulated", "downregulated")) %>%
  dplyr::count(celltype, up_down)
df_nr$celltype <- factor(df_nr$celltype, levels = celltypes)
ggplot(df_nr, aes(x=celltype, y=n)) +
  geom_bar(stat = "identity", fill="darkblue") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab("Cell type") +
  ylab("Number of DE genes (FDR 0.1)") +
  facet_wrap(~up_down, nrow=2) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_1_NumberDEgenes_updown_FDR0.1_", test, "_pcFromCorrectedDataAllCellTypes.pdf"),
       width = 10,
       height = 8)


## number of cells per celltype
cells <- read.csv(paste0(base_workspace, "scanpy_adata/adata_labelTransfer_celltypes_samplesFilt_AnnaAnnotation_obs.csv"))
cells_count <- cells %>%
  dplyr::count(ctAnna_r1) 
cells_count_plot <- left_join(cells_count, df_nr_all, by=c("ctAnna_r1"="celltype"))
ggplot(cells_count_plot, aes(x=n.y, y=n.x, label=ctAnna_r1)) +
  geom_point() +
  geom_text(hjust=-0.15, vjust=-0.15, 
            size=4, color="black",
            check_overlap = TRUE) + 
  xlim(-5,500) +
  xlab("Number of DE genes (FDR 0.1)") +
  ylab("Number of cells") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(filename = paste0(basedir, "plots/04_1_NumberDEgenes_NumberCells_FDR0.1_", test, "_pcFromCorrectedDataAllCellTypes.pdf"),
       width = 10,
       height = 8)


## average fold change for each cell type (up and downregulated)
d_ct <- data.frame(celltype = character(),
                   up_down = character(),
                   mean_fc = numeric())
for (c in celltypes){
  d <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/",
                      c, "_~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
  ) %>%
    dplyr::mutate(up_down = ifelse(log2FoldChange > 0, "upregulated", "downregulated")) %>%
    dplyr::group_by(up_down) %>%
    dplyr::summarise(Mean = mean(log2FoldChange, na.rm=TRUE))
  d_ct <- rbind(d_ct, list("celltype" = c, "up_down" = "downregulated", "mean_fc" = d[d$up_down == "downregulated",]$Mean))
  d_ct <- rbind(d_ct, list("celltype" = c, "up_down" = "upregulated", "mean_fc" = d[d$up_down == "upregulated",]$Mean))
}
d_ct <- left_join(d_ct, df_nr, by=c("celltype","up_down"))
d_ct$n[is.na(d_ct$n)] <- 0
d_ct <- left_join(d_ct, cells_count, by= c("celltype"="ctAnna_r1"))

d_ct$ct_label <- d_ct$celltype
d_ct$ct_label[!d_ct$celltype %in% c("Exc_L2-3", "Exc_L4-6_1", "Exc_L4-6_2",
                                "OPC", "Oligodendrocyte", "In_VIP")] <- NA

d_ct$up_down <- factor(d_ct$up_down, levels = c("upregulated", "downregulated"))

ggplot(d_ct, aes(x=n.x, y=mean_fc, size=n.y, label=ct_label)) +
  geom_point() +
  geom_text(hjust=-0.15, vjust=-0.15, 
            size=4, color="black") + 
  xlim(0,350) +
  xlab("Number of DE genes (FDR 0.1)") +
  ylab("Mean log2(FoldChange)") +
  scale_size_continuous(name="Cluster size") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  ) +
  facet_wrap(~up_down, nrow=2, scales="free")
ggsave(filename = paste0(basedir, "plots/04_1_NumberDEgenesMeanFC_updown_FDR0.1_", test, "_pcFromCorrectedDataAllCellTypes.pdf"),
       width = 8,
       height = 6)


# ## Plots of SLIT2 expression level
# ens_slit2 <- "ENSG00000145147"
# 
# # Read count and metadata
# counts <- readRDS(paste0(basedir, "data/counts.rds"))
# metadata <- readRDS(paste0(basedir, "data/metadata.rds"))
# metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))
# 
# # Split count data to list with entry per celltype
# # Not every cluster is present in all samples; create a vector that represents how to split samples
# splitf <- sapply(stringr::str_split(rownames(counts), 
#                                     pattern = "__",  
#                                     n = 2), 
#                  `[`, 1)
# # Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
# counts_split <- split.data.frame(counts, 
#                                  factor(splitf)) %>%
#   lapply(function(u) 
#     set_colnames(t(u), 
#                  str_split(rownames(u), '__', simplify = TRUE)[,2]))
# 
# slit2 <- lapply(counts_split, function(x) x[ens_slit2,])
# slit2_df <- bind_rows(slit2, .id="celltype") 
# slit2_df_long <- pivot_longer(slit2_df, 
#                               cols = starts_with("SU"),
#                               names_to = "Sample.ID",
#                               values_to = "Reads")
# slit2_df_long <- left_join(slit2_df_long, metadata[c("sample_id", "Sex", "Status")] %>% distinct(), 
#                            by=c("Sample.ID"="sample_id"))
# 
# # plot expression levels across all cell types
# ggplot(slit2_df_long, aes(x=celltype, y=Reads, color=Sex)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(aes(shape=Status),
#              position = position_jitterdodge())
# 
# # plot expression level only in Exc_L4-6_1
# slit2_df_ct <- slit2_df_long %>%
#   filter(celltype == "Exc_L4-6_1")
# ggplot(slit2_df_ct, aes(x=Sex, y=Reads, fill=Status)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge())
# 
# 
# ## Plot fold change of gene Slit2 across all cell types
# d_list <- list()
# for (c in celltypes){
#   d <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/",
#                        c, "_~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch_LRT_pcFromCorrectedDataAllCellTypes.csv")
#   ) %>%
#     dplyr::filter(gene == ens_slit2)
#   d_list[[c]] <- d
# }
# d_df <- bind_rows(d_list, .id = "celltype")
# d_df$celltype <- factor(d_df$celltype, levels = c("Astro_FB", d_df$celltype))
# 
# ggplot(d_df, aes(x=celltype, y=log2FoldChange, size=-log10(padj))) +
#   geom_point() +
#   scale_x_discrete(name = "Cell type",
#                    drop=FALSE,
#                    guide = guide_axis(angle = 45)) +
#   scale_size_continuous(name = "-log10(FDR)") +
#   theme_light() +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 12),
#     strip.text.x = element_text(size = 12)
#   )
# ggsave(filename = paste0(basedir, "plots/04_SLIT2_FC_", test, ".pdf"),
#        width = 8,
#        height = 6)
# 
# 
# # plot RIT2 DE levels across all cell types
# ens_rit2 <- "ENSG00000152214"
# 
# ## Plot fold change of gene Slit2 across all cell types
# d_list <- list()
# for (c in celltypes){
#   d <- read.csv(file.path(basedir, "tables/DESeq2_LRT/",
#                           paste0(c, "_~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch_LRT_pcFromCorrectedDataAllCellTypes.csv"))
#   ) %>%
#     dplyr::filter(gene == ens_rit2)
#   d_list[[c]] <- d
# }
# d_df <- bind_rows(d_list, .id = "celltype")
# d_df$celltype <- factor(d_df$celltype, levels = celltypes)
# 
# ggplot(d_df, aes(x=celltype, y=log2FoldChange, size=-log10(padj))) +
#   geom_point() +
#   scale_x_discrete(name = "Cell type",
#                    drop=FALSE,
#                    guide = guide_axis(angle = 45)) +
#   scale_size_continuous(name = "-log10(FDR)") +
#   theme_light() +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 12),
#     strip.text.x = element_text(size = 12)
#   )
# ggsave(filename = paste0(basedir, "plots/04_RIT2_FC.pdf"),
#        width = 8,
#        height = 6)

