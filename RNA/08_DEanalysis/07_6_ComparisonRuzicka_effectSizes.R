#######################################
### Comparison of effect sizes between our study and Ruzicky SCZ Study ###
#######################################

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(readxl)
library(biomaRt)

## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')
plotdir <- paste0(basedir, "plots/")

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT"

# function to read all sheets of an excel file into list
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# biomart required for mapping of IDs
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

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


# 2.2 Read genes from case control SCZ Ruzicka study
file <-
  paste0(
    basedir,
    "tables/external/Ruzicka_SCZ_Preprint_TableS4.xlsx"
  )
file_data <- read_excel_allsheets(file)


# 3. Correlation of effect sizes
corr_effSize <- function(df_cd, df_scz){
  
  # set gene_symbols as rownames
  rownames(df_cd) <- df_cd$X
  rownames(df_scz) <- df_scz$gene
  
  # subset dataframes to common list of genes
  common_genes <- intersect(rownames(df_cd), rownames(df_scz))
  cor_es <- cor.test(df_cd[common_genes,]$log2FoldChange, 
                     df_scz[common_genes,]$Meta_logFC)
  
  return(c("corr_estimate" = cor_es$estimate,
         "pvalue" = cor_es$p.value))
}

df <- data.frame("ct_cd" = character(),
                 "ct_scz" = character(),
                 "corr_estimate" = numeric(),
                 "pvalue" = numeric())
for (ct_cd in names(all_genes)){
  for(ct_scz in names(file_data)) {
    
    test_comb <- corr_effSize(all_genes[[ct_cd]], file_data[[ct_scz]])
    
    # new_line <- c(ct_cd, ct_scz,
    #               test_comb["corr_estimate.cor"], test_comb["pvalue"])
    # names(new_line) <- c("ct_cd", "ct_scz", "corr_estimate", "pvalue")
    # df <- rbind(df, new_line)
    df <- df %>% add_row("ct_cd" = ct_cd, "ct_scz" = ct_scz,
                         "corr_estimate" = test_comb["corr_estimate.cor"], 
                         "pvalue" = test_comb["pvalue"])
  }
}



# 4 Plot Overlap of cell type hits
ggplot(df, aes(x=ct_cd, y=ct_scz, fill=corr_estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_discrete(drop=FALSE,
                   guide = guide_axis(angle = 45)) +
  xlab("Celltype") +
  ylab("Celltype Ruzicka et al.") +
  theme_light() +
  theme(
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(
  filename = paste0(plotdir, "/07_5_ComparisonHits_ctLevel_Ruzicka/",
                    "/07_5_Barplot_DE_ctLevel_Ruzicka_", test, ".pdf"),
  width = 12,
  height = 8
)  

