## Extract pseudobulk gene score matrices per chromosome

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(RColorBrewer)
library(aod)
library(data.table)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
DEA_dir <- paste0(basedir, "scripts/ATAC/02_DifferentialAnalysis/")

# get argument for chromosome
args = commandArgs(TRUE)
chrom = args[1]
#chrom = "chr1"
#chrom = "chr2"
#chrom = "chr8"

# 1. Read data

projPM <- loadArchRProject(path = "../ArchRSubset_FinalFiltering/")

# 2. Get gene score matrix
gsm <- getMatrixFromProject(
  ArchRProj = projPM,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = c(chrom),
  threads = 1
)


# 3. Generate Pseudo Gene Score matrix per celltype/sample combinations

# add column with pseudobulk sample per cell
gsm@colData$pseudo_group <- paste0(gsm@colData$Sample, "#", gsm@colData$celltypes_refined)

# add pseudobulk ID as colnames
mat_genescore <- gsm@assays@data$GeneScoreMatrix
colnames(mat_genescore) <- gsm@colData$pseudo_group

# sum up gene scores per pseudo group
mat_genescore <- t(mat_genescore)
mat_gsm <- rowsum(mat_genescore, row.names(mat_genescore))

# assign gene names as rownames
colnames(mat_gsm) <- rowData(gsm)$name

# number of cells per pseudosample for normalization
ncells <- table(row.names(mat_genescore))
# normalize pseudosample gene scores
mat_gsm_norm <- sweep(mat_gsm, 1, ncells, "/")
mat_gsm_norm[1:10,1:10]


# 4. Write Pseudobulk gene score matrix to file
fwrite(mat_gsm, file=paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/",
    "Pseudobulk_GeneScoreMatrix_", chrom, ".csv"))
fwrite(mat_gsm_norm, file=paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/",
    "Pseudobulk_GeneScoreMatrix_Normalized_", chrom, ".csv"))

fwrite(data.frame("Pseudosamples" = rownames(mat_gsm)),
      file = paste0(DEA_dir, "tables/Pseudobulk_GeneScoreMatrices/",
                  "Pseudosamples.csv"))