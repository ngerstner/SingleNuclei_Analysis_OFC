## Check the mitochondrial read percentage per cell for QC
## by merging it into ArchR object from singlecell.csv (10X output)

# run with conda environment sc-atac

# 0. load packages
library(ArchR)

# 1. Set reference genome build --> reference used for cellranger is hg38
addArchRThreads(20)
addArchRGenome("hg38")

# 2. Function to merge mitochondrial reads from "singlecell.csv" to ArchR object
#' Get Mitochondrial Reads from singlecell.csv file (10X output) for quality control
#' 
#' @param csvFiles A character vector of names from 10x CSV files to get information regarading mitochondrial reads 
#' @param sampleNames A character vector containing the sample names to be associated with each individual entry in `csvFiles`.
#' @export
getMitoReads10X <- function(
  csvFiles = NULL, 
  sampleNames = NULL
  ){

  if(length(sampleNames) != length(csvFiles)){
    stop("csvFiles and sampleNames must exist!")
  }

  if(!all(file.exists(csvFiles))){
    stop("Not All csvFiles exists!")
  }

  dfList <- lapply(seq_along(csvFiles), function(x){
    df <- suppressPackageStartupMessages(
		suppressMessages(suppressWarnings(readr::read_csv(csvFiles[x])))) 
    if("mitochondrial" %ni% colnames(df)){
      stop("mitochondrial not in colnames of 10x singlecell.csv file! Are you sure input is correct?")
    }
    df <- dplyr::filter(df,
		total > 0)
    df[c("barcode", "total", "mitochondrial")]
  })
  names(dfList) <- sampleNames
  dfConcat <- dplyr::bind_rows(dfList, .id = "sample")
  
  dfConcat
}

# 3. generate list with 10X files
samples <- read.csv("/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/phenotype_clean.csv", header=TRUE)

mitoFileList <- list()
for (mb in c(1,2)) {
	for (i in which(samples$Main.Batch == mb)){
		if (mb == 1){
			basepath <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20201104_Anna_Froehlich_10X_ATACseq/"
		} else {
			basepath <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20210316_Anna_Froehlich_10X_ATACseq/"
		}

		mitoFile <- paste0(basepath,"01_cellranger/SU",samples$SU.Number[i],"_PFC_ATAC/outs/singlecell.csv")
		mitoFileList[paste0("SU",samples$SU.Number[i],"_PFC_ATAC")] <- mitoFile
	}
}
mitoFileList <- unlist(mitoFileList)

# 4. Get number of mitochrondrial reads according to 10x CellRanger
mr <- getMitoReads10X(csvFiles = mitoFileList, sampleNames = names(mitoFileList))
mr <- dplyr::mutate(mr,
	cellID = paste(sample, barcode, sep="#"))

# 5. Merge mitochondrial counts with ArchR project

projPM <- loadArchRProject(path = "ArchROutput/")
cellData <- as.data.frame(projPM@cellColData)
cellData$ID <- rownames(cellData)
mr <- as.data.frame(mr)

df_joined <- dplyr::left_join(cellData, mr, by=c("ID"="cellID"))
df_joined <- dplyr::mutate(df_joined,
			mito_perc = mitochondrial/total)

# 6. Add info to ArchR project and save
projPM <- addCellColData(
  ArchRProj = projPM,
  data = df_joined$mito_perc,
  name = "mito_perc",
  cells = df_joined$ID
)

saveArchRProject(projPM)
