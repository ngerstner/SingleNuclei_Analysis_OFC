## Plot DE hits between PRS extreme groups

# 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)


## 1. Setup

# Define pathes

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# PRS trait
traits <- c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021")

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"

# subsampling to median number of samples (14) in extreme groups
# subsample_mode <- "subsampled"
subsample_mode <- "notSubsampled"

# depending on subsampling mode (yes or no), copy respective files to actual folder
# --> concerns only crossDisorder as it is the only one that gets subsampled
prs_trait <- "crossDisorder2019"
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
            prs_trait, ".csv")
if (subsample_mode == "subsampled") {
  files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext/crossDisorder2019_subsampledToMedianNumberOfSamples"),
                     pattern = p,
                     full.names = TRUE)
  file.copy(files, paste0(basedir, "tables/DESeq2_", test, "/PRSext"),
            overwrite = TRUE)
} else {
  files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext/crossDisorder2019_notSubsampledToMedianNumberOfSamples"),
                      pattern = p,
                      full.names = TRUE)
  file.copy(files, paste0(basedir, "tables/DESeq2_", test, "/PRSext"),
            overwrite = TRUE)
}


## 2. Read in DE genes
# 2.1 Find genes
de_ct <- list()
celltypes <- c()
for (prs_trait in traits){
  p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
              prs_trait, ".csv")
  files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext"),
                      pattern = p,
                      full.names = FALSE)
  # read in the diff genes of different cell types
  de_ct[[prs_trait]] <- list()
  for (f in files){
    de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/PRSext/",f))
    ct <- sub("(\\w*)_~.*","\\1",f)
    de_ct[[prs_trait]][[ct]] <- de
  }
  celltypes <- c(celltypes, names(de_ct[[prs_trait]]))
  de_ct[[prs_trait]] <- de_ct[[prs_trait]][sapply(de_ct[[prs_trait]], 
                                                  function(x) dim(x)[1]) > 0] # remove empty dataframes from list
}
celltypes <- unique(celltypes)

## 2.2. List to dataframe
df_onelevel <- lapply(de_ct, function(x) bind_rows(x, .id="celltype"))
df <- bind_rows(df_onelevel, .id = "prs_trait")
df$celltype <- factor(df$celltype, levels = celltypes)
df$prs_trait <- factor(df$prs_trait, levels = traits)
write.csv(df, paste0(basedir, "tables/DESeq2_", test, "/PRSext/",
                     "04_1_DEhits_allCellTypes.csv"))

## 3. Read in Status DE genes
# find genes
#p <- paste0("X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_ct <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  #de_ct[[ct]] <- de[de$pvalue <= 0.1,]
  de_ct[[ct]] <- de
}
celltypes <- names(de_ct)
de_ct <- de_ct[sapply(de_ct, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# list to dataframe
df_status <- bind_rows(de_ct, .id = "celltype")
df_status$celltype <- factor(df_status$celltype, levels = celltypes)

## 3. Merge PRS DE genes with status DE genes

df_merged <- df %>%
  left_join(df_status, by=c("celltype", "gene")) %>%
  mutate(DEstatus = ifelse(is.na(baseMean.y), FALSE, TRUE))

## 4. Plot DE genes

df_DEstatus <- filter(df_merged, DEstatus == TRUE)
prs_names <- c(
  `crossDisorder2019` = "Cross-Disorder",
  `SCZ2022` = "Schizophrenia",
  `height` = "Height",
  `MDD` = "Major Depressive Disorder",
  `BIP2021` = "Bipolar Disorder"
)
ggplot(df_merged, aes(x=celltype, y=log2FoldChange.x, 
               size = -log10(padj.x))) +
  geom_point(aes(color=prs_trait)) +
  # geom_point(data=df_DEstatus, aes(x=celltype, y=log2FoldChange.x, 
  #                             size = -log10(padj.x)), color = "grey") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  scale_size_continuous(name = "-log10(FDR)") +
  scale_color_manual(name = "GWAS Study",
                     breaks = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021"),
                     labels = c("Cross-Disorder", "Schizophrenia", "Height", "Major Depressive Disorder", "Bipolar Disorder"),
                     values = c("#473335", "#548687", "#CFCFA2", "#B0413E", "#FCAA67")) +
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
ggsave(
  filename = paste0(
    basedir,
    "plots/04_1_DEgenes_FDR0.1_", test, "_PRSadapatednoise_PRSext_propensityScore_allTraits_color_",
    subsample_mode, ".pdf"
  ),
  width = 8,
  height = 6
)


ggplot(df, aes(x=celltype, y=log2FoldChange, 
               size = -log10(padj))) +
  geom_point(aes(color=prs_trait)) +
  geom_point(data=df_DEstatus, aes(x=celltype, y=log2FoldChange.x, 
                                   size = -log10(padj.x)), color = "grey") +
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
    basedir,
    "plots/04_1_DEgenes_FDR0.1_", test, "_PRSadapatednoise_PRSext_propensityScore_allTraits_facet_",
    subsample_mode,
    ".pdf"
  ),
  width = 12,
  height = 8
)


## barplot number of DE genes
df_nr_all <- df %>%
  dplyr::count(celltype, prs_trait)
ggplot(df_nr_all, aes(x=celltype, y=n, fill=prs_trait)) +
  geom_col(position=position_dodge2(preserve = "single")) +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 90)) +
  scale_fill_manual(name = "GWAS Study",
                     breaks = c("crossDisorder2019", "SCZ2022", "height", "MDD", "BIP2021"),
                     labels = c("Cross-Disorder", "Schizophrenia", "Height", "Major Depressive Disorder", "Bipolar Disorder"),
                     values = c("#473335", "#548687", "#CFCFA2", "#B0413E", "#FCAA67"),
                    guide = "none") +
  facet_wrap(~prs_trait, labeller = as_labeller(prs_names)) +#, scales = "free_y") +
  xlab("Cell type") +
  ylab("Number of DE genes (FDR 0.1)") +
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
    basedir,
    "plots/04_1_NumberDEgenes_FDR0.1_", test, "_PRSadaptednoise_PRSext_propensityScore_allTraits_facet_",
    subsample_mode,
    ".pdf"
  ),
  width = 12,
  height = 8
)

