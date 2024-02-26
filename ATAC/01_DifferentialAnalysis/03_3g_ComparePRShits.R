## Comparison of PRS hits between PRS traits based
## on the overlap of genes

# run with conda environment DESeq2

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
    res_list[[ct]] <- res[res$FDR <= 0.1,]
}
res_list <- res_list[sapply(res_list, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# 2.2 Comparison trait
comp_list <- list()
for (ct in celltypes){
    p <- paste0("GeneScore_WaldTest", outlier, norm, "_", ct, "_", prs_trait, ".csv")
    # read in results of differential gene score matrix analysis
    res <- fread(paste0(DEA_dir, "tables/GeneScoreMatrix/PRSext/", p))
    comp_list[[ct]] <- res[res$FDR <= 0.1,]
}
comp_list <- comp_list[sapply(comp_list, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# 3. Compare hits per celltype
df_plot <- data.frame(celltype = character(),
                      overlap = numeric(),
                      comp_only = numeric(),
                      height_only = numeric())
for (ct in celltypes){
  
  # number of genes present in both w and wo status covariate model 
  overlap_w_wo <- length(intersect(res_list[[ct]]$gene, comp_list[[ct]]$gene))
  only_wo <-  length(res_list[[ct]]$gene) - overlap_w_wo
  only_w <-  length(comp_list[[ct]]$gene) - overlap_w_wo
  
  # add numbers to dataframe for plotting
  df_plot <- rbind(df_plot,
                   list("celltype" = ct,
                        "overlap" = overlap_w_wo, "comp_only" = only_w,
                        "height_only" = only_wo))
}


# 4. Plot Comparsion of hits for model w and wo Status covariate
df_plot <- df_plot %>%
  pivot_longer(cols = overlap:height_only, names_to = "category") %>%
  mutate_at(c("celltype", "category"), factor) %>%
  mutate(category = fct_relevel(category, "height_only", "overlap", "comp_only"))

ggplot(df_plot, aes(x=celltype, y=value, fill=category)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 90)) +
  scale_fill_viridis(discrete = TRUE, name = "",
                     labels = c(paste0(ref_trait, " - DA risk genes"), "Overlap",
                                paste0(prs_trait, " - DA risk genes"))) +
  xlab("Cell type") +
  ylab("Number of sig. DA genes") +
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
  filename = paste0(
    DEA_dir,
    "plots/03_3g_ComparePRShits/03_3g_ComparePRShits_",ref_trait,"-", prs_trait, ".pdf"
  ),
  width = 12,
  height = 8
)  



