## Propensity score matching for PRS extreme groups

## 0. Load Packages

# load edgeR which loads limma as dependency
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(MatchIt)
library(optmatch)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2/')
prsdir <- '~/Documents/PostmortemBrain/genotype/brainbank_PRS/'

# source file with functions needed
source(paste0(basedir,"00_functions.R"))

# celltype - not really necessary but code was established like this
cluster <- "Exc_L2-3"

# PRS trait
# prs_trait <- "crossDisorder2019"
# prs_trait <- "MDD"
# prs_trait <- "SCZ2022"
# prs_trait <- "BIP2021"
prs_trait <- "height"


# cutoff --> which proportion of sample needs to be expressed for gene to be kept?
cutoff_sample <- 0.75


# # numbers of samples per extreme group across PRS traits
# # --> subsample to median of these numbers 
# n_extGroups <- c("BIP2021" = 14, "crossDisorder2019" = 17,
#                       "height" = 14, "MDD" = 11, "SCZ2022" = 13)
# med_extGroups <- median(n_extGroups)


# 2. Preprocessing of data

# Read count and metadata
metadata <- readRDS(paste0(basedir, "data/metadata.rds"))
metadata$Mode.of.death <- as.factor(str_replace_all(metadata$Mode.of.death, " ", ""))


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


# read PRS from specified trait 
prs <- read.table(paste0(prsdir,prs_trait,"_plink_out.sscore")) %>%
  mutate(PRS = scale(V5)) %>%
  mutate(sample_id = paste0(V2,"_PFC_RNA"))


# join metadata with PRS
cluster_metadata <- cluster_metadata %>%
  left_join(prs[c("sample_id","PRS")], by="sample_id")


# 3. Get groups with high and low polygenic risk using propensity score matching

# number of samples for 15% extremes
n_ext <- round(0.15*nrow(cluster_metadata))

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
                     # caliper = 1,
                     replace = FALSE)
summary(match_obj)

plot(match_obj, type = "jitter", interactive = FALSE)
plot(summary(match_obj), abs = FALSE)

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


# 4. Plot characteristics of high and low group
# Age
p_age <- ggplot(matched_data, aes(x = ExtGroup, y=Age, fill = ExtGroup)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave(p_age, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                prs_trait, "_Age.pdf"),
       width = 10,
       height = 8)

# Sex
p_sex <- ggplot(matched_data, aes(x = ExtGroup, fill = Sex)) +
  geom_bar(stat="count") + 
  xlab("PRS group") +
  scale_fill_brewer(palette="Set1")+
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave(p_sex, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                prs_trait, "_Sex.pdf"),
       width = 10,
       height = 8)

# Brain.pH
p_ph <- ggplot(matched_data, aes(x = ExtGroup, y=Brain.pH, fill = ExtGroup)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave(p_ph, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                               prs_trait, "_Brain.pH.pdf"),
       width = 10,
       height = 8)


# PMI
p_pmi <- ggplot(matched_data, aes(x = ExtGroup, y=PMI, fill = ExtGroup)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave(p_pmi, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                prs_trait, "_PMI.pdf"),
       width = 10,
       height = 8)

# RIN
p_rin <- ggplot(matched_data, aes(x = ExtGroup, y=RIN, fill = ExtGroup)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave(p_rin, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                prs_trait, "_RIN.pdf"),
       width = 10,
       height = 8)

# X6.Batch
p_batch <- ggplot(matched_data, aes(x = ExtGroup, fill = X6.Batch)) +
  geom_bar(stat="count") + 
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave(p_batch, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                  prs_trait, "_X6.Batch.pdf"),
       width = 10,
       height = 8)

# Case/Control status
p_cc <- ggplot(matched_data, aes(x = ExtGroup, fill = Status)) +
  geom_bar(stat="count") + 
  scale_fill_manual(breaks = c(0,1),
                    labels = c("Control", "Case"),
                    values = c("#EFEE68", "#48A4B6")) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave(p_batch, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                  prs_trait, "_Status.pdf"),
       width = 10,
       height = 8)

# Disease classification
p_class <- ggplot(matched_data, aes(x = ExtGroup, fill = Classification)) +
  geom_bar(stat="count") + 
  scale_fill_manual(breaks = c("BipolarDisorder", "Control", "MajorDepression",
                               "SchizoaffectiveDisorder", "Schizophrenia"),
                    # labels = c("Control", "Case"),
                    values = c("#773D44","#806A8E","#48A4B6","#65D498","#EFEE68")) +
  xlab("PRS group") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave(p_class, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                                  prs_trait, "_Classification.pdf"),
       width = 10,
       height = 8)

# save all plots in one file
ggarrange(p_age, p_ph, p_sex, p_cc, p_rin, p_pmi, p_batch,p_class,
          ncol = 4, nrow = 2)
ggsave(filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_",
                         prs_trait, "_all.pdf"),
       width = 16,
       height = 12)


# histogram of PRS
p_hist <- ggplot() +
  geom_histogram(data=cluster_metadata, mapping=aes(x=PRS)) +
  geom_histogram(data=matched_data, mapping=aes(x=PRS, fill=ExtGroup)) +
  ggtitle(paste0("High: ", table(matched_data$ExtGroup)[1], ", Low: ", table(matched_data$ExtGroup)[2])) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave(p_hist, filename = paste0(basedir, "plots/01_PRS_extremeGroups/propensityScore/01_PRS_extremeGroups_propensityScore_", 
                                 prs_trait, "_histogram.pdf"),
       width = 10,
       height = 8)


# write samples to file
write.csv(matched_data, file = paste0(basedir, "tables/PRSext/samples_propensityScore/",
                                      prs_trait,"_samples_propensityScore.csv"))
