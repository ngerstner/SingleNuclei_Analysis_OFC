## Differential expression analysis between PRS extreme groups
## with DESeq2 on pseudobulk data

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
library(MatchIt)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')
prsdir <- '~/Documents/PostmortemBrain/genotype/brainbank_PRS/'

# source file with functions needed
source(paste0(basedir,"00_functions.R"))

# celltype of interest
args = commandArgs(TRUE)
cluster = args[1]
# cluster <- "Exc_L2-3"
# cluster <- "Oligodendrocyte"
# cluster <- "OPC"
# cluster <- "Astro_FB"
# cluster <- "Astro_PP"
# cluster <- "Exc_L3-5"
# cluster <- "Exc_L4-6_1"
# cluster <- "In_LAMP5"
# cluster <- "In_VIP"

# PRS trait
# prs_trait <- "crossDisorder2019"
# prs_trait <- "MDD"
# prs_trait <- "SCZ2022"
# prs_trait <- "BIP2021"
prs_trait <- "height"

# test: Likelihood Ratio Test or Wald Test
test <- "Wald"
# test <- "LRT"

# plots for top hits
plot_th <- TRUE
# plot_th <- FALSE

# cutoff --> which proportion of sample needs to be expressed for gene to be kept?
cutoff_sample <- 0.75

# numbers of samples per extreme group across PRS traits
# --> subsample to median of these numbers
n_extGroups <- c("BIP2021" = 14, "crossDisorder2019" = 17,
                 "height" = 14, "MDD" = 11, "SCZ2022" = 13)
med_extGroups <- median(n_extGroups)

# des <- "~ PRSext + Sex + Age + Brain.pH + RIN + PMI + X6.Batch"
# des_null <- "~ Sex + Age + Brain.pH + RIN + PMI + X6.Batch"
# des_varPart <- "~ (1|PRSext) + (1|Sex) + Age + Brain.pH + RIN + PMI + (1|X6.Batch)"

des <- "~ PRSext + Sex + Age + Brain.pH + RIN + PMI"
des_null <- "~ Sex + Age + Brain.pH + RIN + PMI"
des_varPart <- "~ (1|PRSext) + (1|Sex) + Age + Brain.pH + RIN + PMI"


# 2. Preprocessing of data

# Read count and metadata
counts <- readRDS(paste0(basedir, "data/counts.rds"))
metadata <- readRDS(paste0(basedir, "data/metadata.rds"))
metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))

# Read list of mitochondrial genes
mt_genes <- read.csv(paste0(base_workspace, "/tables/mito_genes.csv"))

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

# read PC1 from inferred from corrected data in all cell types
pc <- read.csv(paste0(basedir,"/tables/PCA_noise/All_PCAnoise_rotated_PRS_",prs_trait,".csv"))

# read PRS from specified trait 
prs <- read.table(paste0(prsdir,prs_trait,"_plink_out.sscore")) %>%
  mutate(PRS = scale(V5)) %>%
  mutate(sample_id = paste0(V2,"_PFC_RNA"))


# join metadata with PC1 from corrected data and PRS
cluster_metadata <- cluster_metadata %>%
  left_join(as.data.frame(pc[c("X", "PC1")]),
            by = c("sample_id" = "X")) %>%
  left_join(prs[c("sample_id","PRS")], by="sample_id")

cluster_metadata <- dplyr::rename(cluster_metadata,
                                  PCnoise = PC1)

# order samples by PRS
cluster_metadata <- arrange(cluster_metadata, desc(PRS))

# get extreme subsets of metadata
metadata_high20 <- head(cluster_metadata, n = 20) %>%
  mutate(ExtGroup = "high")
metadata_low20 <- tail(cluster_metadata, n = 20) %>%
  mutate(ExtGroup = "low")
metadata_ext <- rbind(metadata_high20, metadata_low20)
metadata_ext$ExtGroup <- as.factor(metadata_ext$ExtGroup)


# match each sample from high group to one of low group
match_obj <- matchit(ExtGroup ~ Age + Sex + Brain.pH + PMI + RIN,
                     data = metadata_ext, method = "nearest", distance ="glm",
                     exact = ~Sex,
                     discard = "both",
                     ratio = 1,
                     replace = FALSE)

# Extract the matched data and save the data into the variable matched_data
matched_data <- match.data(match_obj)

# # if more than median number of samples per extreme group (here 14) --> subsample to 14 per group
# if (nrow(matched_data) > 2*med_extGroups){
#   n_perGroup <- nrow(matched_data)/2
#   n_removePerGroup <- n_perGroup - med_extGroups
# 
#   # select n_removePerGroup matched pairs and remove them
#   groups <- unique(matched_data$subclass)
#   set.seed(8)
#   indices_remove <- sample(groups, n_removePerGroup)
#   matched_data <- filter(matched_data, !(subclass %in% indices_remove))
# }

#Extract the matched data and save the data into the variable matched_data
cluster_metadata_filt <- matched_data %>%
  dplyr::rename(PRSext = ExtGroup)

# subset cluster_counts
cluster_counts_filt <- cluster_counts[,cluster_metadata_filt$sample_id]

# DESeq object
dds <- DESeqDataSetFromMatrix(cluster_counts_filt, 
                              colData = cluster_metadata_filt, 
                              design = as.formula(des))
dds

# Remove mitochondrial genes
dds <- dds[!rownames(dds) %in% mt_genes$gene_ids,]

# Keep only genes expressed in more than specific percentage of samples
dds <- dds[rowSums(counts(dds) >= 10) >= cutoff_sample*ncol(dds), ]
dds

# Remove outliers
dds <- outlier_removal(dds, paste0(basedir, "plots/outlier_removal/",prs_trait, 
                                   "/03_1_PRSext_propensityScore_",
                                   prs_trait,"_",cluster))

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

# get results
# low is reference level --> foldchange refers to change from low to high
res <- results(dds,
               contrast = c("PRSext", "high", "low"),
               # name = "PRSext_low_vs_high",
               alpha = 0.1)
res <- lfcShrink(dds,
                 contrast = c("PRSext", "high", "low"),
                 # coef = "PRSext_low_vs_high",
                 res=res,
                 type = "normal")

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Set thresholds
padj_cutoff <- 0.1
# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
# join with ENSEMBL IDs
genes <- read.csv(paste0(base_workspace, "tables/geneIDs_pseudocounts.csv"))
sig_res <- left_join(sig_res, genes, by=c("gene" = "gene_ids"))
res_tbl <- left_join(res_tbl, genes, by=c("gene" = "gene_ids"))

# Check significant genes output
sig_res

# Write significant hits to file
write.csv(sig_res,
          paste0(basedir, "tables/DESeq2_", test, "/PRSext/", cluster, "_", des, "_sig", padj_cutoff, 
                 "_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", prs_trait, ".csv"))

# Write all genes to file
write.csv(res_tbl,
          paste0(basedir, "tables/DESeq2_", test, "/PRSext/", cluster, "_", des,
                 "_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", prs_trait, ".csv"))


# metadata DESeq analysis
meta_deseq <- data.frame("samples_total" = nrow(colData(dds)),
                         "genes_tested" = nrow(counts(dds))
                         )
# Write all genes to file
write.csv(meta_deseq,
          paste0(basedir, "tables/DESeq2_", test, "/PRSext/", cluster, "_", des,
                 "_", test, "_pcFromCorrectedDataAllCellTypes_PRSadaptedNoise_PRSext_propensityScore_", prs_trait, "_metadata.csv"))




# 5. Plot expression levels of top hits

# TODO: run again so that outlier points are not plotted from boxplot
# they are plotted by jitter again

# function to plot tophits
tophit_boxplots <- function(g, gene_name, data = c("raw", "norm", "corr")){
  
  # subset counts to gene
  if(data == "raw"){
    counts_tophit <- as.data.frame(counts(dds)[g,])
    colnames(counts_tophit) <- c(g)
  } else if (data == "norm") {
    counts_tophit <- as.data.frame(assay(vsd)[g,])
    colnames(counts_tophit) <- c(g)
  } else {
    y3 <- removeBatchEffect(y, batch = d$samples$Sex, # batch2 = d$samples$X6.Batch,
                        covariates = colData(dds)[,c("Age", "RIN", "Brain.pH", "PMI", "PCnoise")])
    counts_tophit <- as.data.frame(y3[g,])
    colnames(counts_tophit) <- c(g)
  }
  counts_tophit <- mutate(counts_tophit, sample = rownames(counts_tophit))
  # counts_tophit <- left_join(counts_tophit, cluster_metadata, by = c("sample" = "sample_id"))
  counts_tophit <- left_join(counts_tophit, as.data.frame(colData(dds)), by = c("sample" = "sample_id"))
  
  # boxplot status
  ggplot(counts_tophit, aes_string(x="PRSext", y=g, fill = "PRSext")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter() +
          scale_fill_brewer(palette="Dark2") +
          ylab("Expression level") +
          ggtitle(paste0(g, " - ", gene_name))
          theme_light() +
          theme(
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 12)
          )
  ggsave(filename = paste0(basedir, "plots/03_1_DESeq2_tophits/", prs_trait, "/",cluster,"_", gene_name, "_", data, ".pdf"),
       width = 8,
       height = 6)
  
}  

if (plot_th) {
  
  # correct for batch effects with limma
  # Create DGEList object
  d0 <- DGEList(counts(dds), samples = as.data.frame(colData(dds)))
  dim(d0)
  
  # Filter genes
  cutoff_counts <- 10
  keep <- rowSums(d0$counts >= cutoff_counts) >= cutoff_sample*ncol(d0$counts)
  d <- d0[keep,] 
  dim(d) # number of genes left
  
  # Normalization factor
  d <- calcNormFactors(d)
  d
  
  # specify model to fitted
  mm <- model.matrix(as.formula(des), data = d$samples)
  # voom transformation
  y <- voom(d, mm, plot = T)
  
  # plot some of the tophits
  if (nrow(sig_res) >= 1){
  tophit_boxplots(sig_res$gene[1], sig_res$X[1], data = "corr")
  tophit_boxplots(sig_res$gene[1], sig_res$X[1], data = "norm")
  tophit_boxplots(sig_res$gene[1], sig_res$X[1], data = "raw")
  }
  
  if (nrow(sig_res) >= 2){
  tophit_boxplots(sig_res$gene[2], sig_res$X[2], data = "corr")
  tophit_boxplots(sig_res$gene[2], sig_res$X[2], data = "norm")
  tophit_boxplots(sig_res$gene[2], sig_res$X[2], data = "raw")
  }
  
  if (nrow(sig_res) >= 3){
  tophit_boxplots(sig_res$gene[3], sig_res$X[3], data = "corr")
  tophit_boxplots(sig_res$gene[3], sig_res$X[3], data = "norm")
  tophit_boxplots(sig_res$gene[3], sig_res$X[3], data = "raw")
  }
  
  
  
  
  
}
