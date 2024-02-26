## Comparison of PRS hits between traits looking
## at the overlap between DE hits

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(viridis)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# reference trait to compare to 
ref_trait <- "height"
# ref_trait <- "crossDisorder2019"
# ref_trait <- "MDD"
# ref_trait <- "SCZ2022"

# PRS trait
# prs_trait <- "crossDisorder2019"
# prs_trait <- "MDD"
# prs_trait <- "SCZ2022"
prs_trait <- "BIP2021"

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"

## 1. Read in DE genes
# 1.1 Height
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
            ref_trait, ".csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_height <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/PRSext/",f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_height[[ct]] <- de
}
de_height <- de_height[sapply(de_height, function(x) dim(x)[1]) > 0] # remove empty dataframes from list

# 2.2 Comparison trait
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", 
            prs_trait, ".csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/PRSext"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_comp <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/PRSext/",f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_comp[[ct]] <- de
}
celltypes <- names(de_comp)
de_comp <- de_comp[sapply(de_comp, function(x) dim(x)[1]) > 0] # remove empty dataframes from list


# 3. Compare hits per celltype
df_plot <- data.frame(celltype = character(),
                      overlap = numeric(),
                      comp_only = numeric(),
                      height_only = numeric())
for (ct in celltypes){
  
  # number of genes present in both w and wo status covariate model 
  overlap_w_wo <- length(intersect(de_height[[ct]]$gene, de_comp[[ct]]$gene))
  only_wo <-  length(de_height[[ct]]$gene) - overlap_w_wo
  only_w <-  length(de_comp[[ct]]$gene) - overlap_w_wo
  
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
                     labels = c(paste0(ref_trait, " - DE risk genes"), "Overlap",
                                paste0(prs_trait, " - DE risk genes"))) +
  xlab("Cell type") +
  ylab("Number of sig. DE genes") +
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
    basedir,
    "plots/04_4_ComparePRShits/04_4_ComparePRShits_",ref_trait,"-", prs_trait, ".pdf"
  ),
  width = 12,
  height = 8
)  



