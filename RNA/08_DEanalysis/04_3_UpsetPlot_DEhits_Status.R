## Generate UpSet Plot to visualize the overlap 
## of DE hits between cell types

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(UpSetR)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"


## 2. Read in DE genes
p <- paste0("sig0.1_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
# read in the diff genes of different cell types
de_ct <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  de_ct[[ct]] <- de$gene
}
celltypes <- names(de_ct)


## 3. Make UpSet Plot 
pdf(file = paste0(basedir, "/plots/04_3_UpsetPlot_DEhits_Status.pdf"), 
    height = 7, width = 10)
#par(mar = c(5, 10, 4, 2) + 0.1)
print(upset(fromList(de_ct), 
            sets = rev(names(de_ct)), nsets = length(de_ct), 
            keep.order = TRUE, order.by = "freq",
            mb.ratio = c(0.4, 0.6),
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8),
            set_size.show = TRUE,
            sets.x.label = "#DE genes in celltype",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()


## 4. Make UpSet Plot without the celltypes with 0 hits
de_ct <- de_ct[ sapply(de_ct, length) > 0 ]
pdf(file = paste0(basedir, "/plots/04_3_UpsetPlot_DEhits_Status_nonEmpty.pdf"), 
    height = 7, width = 10)
#par(mar = c(5, 10, 4, 2) + 0.1)
print(upset(fromList(de_ct), 
            sets = rev(names(de_ct)), nsets = length(de_ct), 
            keep.order = TRUE, order.by = "freq",
            mb.ratio = c(0.4, 0.6),
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8),
            set_size.show = TRUE,
            sets.x.label = "#DE genes in celltype",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()
