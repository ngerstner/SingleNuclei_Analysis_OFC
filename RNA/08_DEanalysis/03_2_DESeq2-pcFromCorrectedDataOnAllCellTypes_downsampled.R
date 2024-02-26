## Differential expression analysis between disease groups
## with DESeq2 on downsampled pseudobulk data

## 0. Load Packages

# load edgeR which loads limma as dependency
library(DESeq2)
library(edgeR)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(variancePartition)
library(PCAtools)
library(lattice)
library(sva)


## 1. Setup

# Define pathes, cell type of interest and model.

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# source file with functions needed
source(paste0(basedir,"00_functions.R"))

# celltype of interest
args = commandArgs(TRUE)
cluster = args[1]
# cluster <- "Exc_L2-3"
# cluster <- "Oligodendrocyte"
# cluster <- "OPC"
# cluster <- "Astro_FB"
# cluster <- "Exc_L3-5"
# cluster <- "Exc_L4-6_1"
# cluster <- "In_LAMP5"

# percentile of downsampling
perc = as.numeric(args[2])
# perc = 75

# test: Likelihood Ratio Test or Wald Test
test <- "Wald"
# test <- "LRT"

# cutoff --> which proportion of sample needs to be expressed for gene to be kept?
cutoff_sample <- 0.75

des <- "~ Status + Sex + Age + Brain.pH + RIN + PMI + X6.Batch"
des_null <- "~ Sex + Age + Brain.pH + RIN + PMI + X6.Batch"
des_varPart <- "~ (1|Status) + (1|Sex) + Age + Brain.pH + RIN + PMI + (1|X6.Batch)"


# 2. Preprocessing of data

# Read count and metadata
counts <- readRDS(paste0(basedir, "data/counts_perc", perc, ".rds"))
metadata <- readRDS(paste0(basedir, "data/metadata_perc", perc, ".rds"))
metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))

# 1a. Split count data to list with entry per celltype
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
# Explore the different components of list
str(counts_split)


# 2. Prepare data
# Scale and center numeric variables with very different range
metadata$cell_count_sample <- scale(metadata$cell_count_sample)
metadata$cell_count_pseudosample <- scale(metadata$cell_count_pseudosample)
# Subset the metadata to only the cluster of interest
cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the B cells
counts <- counts_split[[cluster]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))


# 3. Create DESeq object

# Create DESeq object that includes also PC1 inferred from the corrected data in all
# cell types as covariate. 
# Filter genes according to cutoff.
# Remove all samples that are more than 3SDs away from mean in PC1 (outliers).
# Run variance stabilizing transformation.

# read PC1 from inferred from corrected data in all cell types
pc <- read.csv(paste0(basedir,"/tables/PCA_noise/All_PCAnoise_perc", perc, "_rotated.csv"))

# join metadata with PC1 from corrected data
cluster_metadata <- left_join(cluster_metadata, 
                              as.data.frame(pc[c("X", "PC1")]), 
                              by=c("sample_id"="X"))
cluster_metadata <- rename(cluster_metadata,
                           PCnoise = PC1)

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = as.formula(des))
dds

# Keep only genes expressed in more than specific percentage of samples
dds <- dds[rowSums(counts(dds) >= 10) >= cutoff_sample*ncol(dds), ]
dds

# Remove outliers
dds <- outlier_removal(dds, paste0(basedir, "plots/outlier_removal/downsampled",
                                 "/03_2_Downsampling",
                                 perc,"_",cluster))

# Normalization of count data
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
hist(assay(vsd))


# 4. Run Differential Expression
  
# Run DESeq2 differential expression analysis
design(dds) <- as.formula(paste0(des,"+PCnoise"))

if (test == "LRT"){
  # a) with LRT
  # reduced model is the model without the variable of interest --> ALSO FOR SVA?
  dds <- DESeq(dds, test = "LRT", reduced = as.formula(paste0(des_null,"+PCnoise")))
} else {
  # b) with Wald test
  dds <- DESeq(dds)
}

# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast for stim vs ctrl
levels(cluster_metadata$Status)[2]
levels(cluster_metadata$Status)[1]

# get results
contrast <- c("Status", levels(cluster_metadata$Status)[2], levels(cluster_metadata$Status)[1])
res <- results(dds, 
               contrast = contrast,
               alpha = 0.1)
res <- lfcShrink(dds, 
                 contrast = contrast,
                 res=res,
                 type = "normal")

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl
write.csv(res_tbl,
          paste0(basedir, "tables/DESeq2_", test, "/downsampled/", cluster, "_", des, 
                 "_", test, "_pcFromCorrectedDataAllCellTypes_perc", perc, ".csv"))


# Set thresholds
padj_cutoff <- 0.1
# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
# join with ENSEMBL IDs
genes <- read.csv(paste0(base_workspace, "tables/geneIDs_pseudocounts.csv"))
sig_res <- left_join(sig_res, genes, by=c("gene" = "gene_ids"))

# Check significant genes output
sig_res

# Write significant hits to file
write.csv(sig_res, 
          paste0(basedir, "tables/DESeq2_", test, "/downsampled/", cluster, "_", des, 
                 "_sig", padj_cutoff, "_", test, "_pcFromCorrectedDataAllCellTypes_perc",perc,".csv"))


# metadata DESeq analysis
meta_deseq <- data.frame("samples_total" = nrow(colData(dds)),
                         "genes_tested" = nrow(counts(dds))
)
# Write all genes to file
write.csv(meta_deseq,
          paste0(basedir, "tables/DESeq2_", test, "/downsampled/", cluster, "_", des,
                 "_", test, "_pcFromCorrectedDataAllCellTypes_perc", perc, "_metadata.csv"))




