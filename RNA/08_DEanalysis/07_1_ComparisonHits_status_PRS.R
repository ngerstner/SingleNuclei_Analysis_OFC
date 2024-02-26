## Compare the disease status DE hits to the PRS extreme group
## DE hits 

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggpubr)

## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')
plotdir <- paste0(basedir, "plots/")

# PRS traits
prs_traits <- c("crossDisorder2019", "MDD", "SCZ2022", "BIP2021", "height")

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"

## 2. Read in DE genes
de_genes <- list()
all_genes <- list()

# 2.1 Read genes from case control
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
  de_genes[[ct]] <- list()
  de_genes[[ct]][["diseaseStatus"]] <- de[de$padj <= 0.1,]
  all_genes[[ct]] <- list()
  all_genes[[ct]][["diseaseStatus"]] <- de
  # print(ct)
  # print(nrow(de[de$padj <= 0.1,]))
}
celltypes <- names(de_genes)

# 2.2 Read genes from PRSext
for (prs_trait in prs_traits){
  p <- paste0("PMI_",test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_",prs_trait, ".csv")
  files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext/"),
                      pattern = p,
                      full.names = FALSE)
  files <- files[!str_detect(files, "metadata")]
  # read in diff genes of different cell types
  for (f in files){
    de <- read.csv(paste0(basedir, "tables/DESeq2_",test, "/PRSext/",f))
    de[is.na(de)] <- 1
    ct <- sub("(\\w*)_~.*","\\1",f)
    de_genes[[ct]][[prs_trait]] <- de[de$padj <= 0.1,]
    all_genes[[ct]][[prs_trait]] <- de
    # print(nrow(de[de$padj <= 0.1,]))
  }
}


# 3. Plots

# 3.1 Overlap of hits per cell type
for (ct in celltypes){
  
  l <- list("diseaseStatus" = de_genes[[ct]][["diseaseStatus"]]$gene,
            "crossDisorder2019" = de_genes[[ct]][["crossDisorder2019"]]$gene,
            "MDD" = de_genes[[ct]][["MDD"]]$gene,
            "SCZ" = de_genes[[ct]][["SCZ"]]$gene,
            "BIP2021" = de_genes[[ct]][["BIP2021"]]$gene)
  l <- l[lapply(l,length)>0]
  
  if (length(l) > 0){
    m <- ComplexHeatmap::make_comb_mat(l)
    print(m)
    pdf(file = paste0(plotdir, "/07_1_ComparisonHits_status_PRS/",
                      "/07_1_",ct,"_DE_PRS_", test, ".pdf"),
        width = 12, height = 10)
    print(ComplexHeatmap::UpSet(m))
    dev.off()
  }
}

# 3.2 Comparison of FC and p-values in other comparison to hit comparison
for (ct in celltypes){
  
  # subset lists to respective celltype
  ct_de <- de_genes[[ct]]
  ct_all <- all_genes[[ct]]
  
  # names of comparisons, e.g. diseaseStatus, crossDisorder2019
  comp_names <- names(ct_all)
  
  # reformat lists to dataframes
  ct_de <- bind_rows(ct_de, .id="comparison")
  ct_all <- bind_rows(ct_all, .id="comparison")
  
  # add FC and p-values of other comparisons to DE table
  ct_all <- ct_all %>%
    select(comparison, gene, log2FoldChange, padj) %>%
    filter(gene %in% ct_de$gene) %>%
    mutate(de_cutoff = (padj<=0.1))
  ct_de <- ct_de %>%
    select(comparison, gene, log2FoldChange, padj) %>%
    left_join(ct_all, by="gene")
  n_genes <- ct_de %>%
    filter(de_cutoff) %>%
    group_by(comparison.x, comparison.y) %>%
    count(de_cutoff) 
  ct_de <- ct_de %>%
    mutate(padj.y = -log10(padj.y)) %>%
    pivot_longer(cols = c("log2FoldChange.y", "padj.y"))
  
  # change comparison.x to factor with all comparisons as levels
  ct_de$comparison.x <- factor(ct_de$comparison.x, levels = comp_names)
  ct_de$comparison.y <- factor(ct_de$comparison.y, levels = comp_names)
  
  # change comparison.x to factor with all comparisons as levels
  n_genes$comparison.x <- factor(n_genes$comparison.x, levels = comp_names)
  n_genes$comparison.y <- factor(n_genes$comparison.y, levels = comp_names)
  # expand dataframe to include all possible pairs of comparisons
  n_genes <- n_genes %>%
    dplyr::select(-de_cutoff) %>%
    as.data.frame() %>%
    tidyr::complete(comparison.x, comparison.y) %>%
    distinct() %>%
    mutate(label = factor("Number of sig. DE genes"))
  
  if(nrow(ct_de) > 0){
  # boxplot showing hits within the ct detected and in the other cell types
  g_box <- ggplot(ct_de) +
    geom_boxplot(aes(x = comparison.x, y = value, fill = comparison.y, col = comparison.y)) +
    facet_wrap(~name, nrow = 2, scales = "free_y",
               labeller = labeller(name = c("log2FoldChange.y"= "log2FoldChange",
                                            "padj.y" = "-log10(FDR)",
                                            "n" = "Number of DE genes"))) +
    scale_x_discrete(drop=FALSE,
                     guide = guide_axis(angle = 45)) +
    scale_color_brewer(palette="Paired") +
    scale_fill_brewer("", palette="Paired") +
    xlab("DE genes were identified in") +
    ylab("") +
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
  # add median line in black
  dat <- ggplot_build(g_box)$data[[1]]
  dat$name <- as.factor(c(rep("log2FoldChange.y", nrow(dat)/2), 
                          rep("padj.y", nrow(dat)/2)))
  g_box <- g_box + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                    y=middle, yend=middle), colour="black", size=0.25,
                       inherit.aes = FALSE) 
  
  # barplot showing numer of DE genes overlapping with other comparison
  g_bar <- ggplot(n_genes, aes(x = comparison.x, y = n, fill = comparison.y)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_x_discrete(drop=FALSE,
                     guide = guide_axis(angle = 45)) +
    scale_color_brewer(palette="Paired") +
    scale_fill_brewer("", palette="Paired") +
    xlab("DE genes were identified in") +
    ylab("") +
    facet_wrap(~label) +
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
  
  # combine boxplot and barplot
  ggarrange(g_bar + 
              theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank() ), 
            g_box,
            ncol = 1,
            align = "v",
            heights = c(1,2.5),
            legend = "bottom",
            common.legend = TRUE)
  ggsave(
    filename = paste0(plotdir, "/07_1_ComparisonHits_status_PRS/",
                      "/07_1_Boxplot_",ct,"_DE_PRS_", test, ".pdf"),
    width = 12,
    height = 8
  )  
  }
}

