## Add GeneIntegrationMatrix to ArchR project
## each cell gets gene expression values from "closest neighbour" in RNA data
## and add GeneIntegrationMatrix to arrow files

# run with conda environment sc-atac

library(ArchR)
library(Seurat)
library(SeuratDisk)

# 1. Convert anndata object to Seurat object (only necesarry once)
#Convert(source = "../../scanpy_adata/adata_labelTransfer_celltypes_samplesFilt_AnnaAnnotation_woObs.h5ad", 
#    dest="h5seurat", 
#    assay = "RNA", 
#    overwrite=TRUE)

# 2. Read Seurat object containing RNAseq data
rna <- LoadH5Seurat("../../scanpy_adata/adata_labelTransfer_celltypes_samplesFilt_AnnaAnnotation_woObs.h5seurat",
    assays="RNA")

# 3. Read metadata/adata.obs as conversion from anndata --> Seurat works only without adata.obs
metadata <- fread("../../scanpy_adata/adata_labelTransfer_celltypes_samplesFilt_AnnaAnnotation_obs.csv")

# 4. Assign metadata dataframe to meta.data slot of object
# IMPORTANT: row.names need to be the same (also order) as column names of counts for the integration
rna@meta.data <- as.data.frame(metadata)
rownames(rna@meta.data) <- rna@meta.data$index
#format(object.size(rna), units="GB")
#[1] "38.9 Gb"

# 5. Load ArchR project with ATACseq data
projPM <- loadArchRProject(path = "ArchROutput/")
#format(object.size(projPM), units="MB")
#[1] "922.5 Mb"

# deactivate parallelization as unstable here
addArchRThreads(1)
# 6. Add Pseudo-scRNA-seq profiles for each ATAC cell
projPM <- addGeneIntegrationMatrix(
    ArchRProj = projPM, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna,
    addToArrow = TRUE,
    force = TRUE,
    groupRNA = "Anna_celltypes",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

# 7. Add imputation weights (describe each cell as a linear combination of many cells based on a MAGIC diffusion matrix)
projPM <- addImputeWeights(projPM)

# 8. Save ArchR project 
saveArchRProject(ArchRProj = projPM)		
