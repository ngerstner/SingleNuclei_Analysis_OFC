## Create arrow files from ATAC fragment files (part of CellRanger output)

# run with conda environment sc-atac

# 0. load packages
library(ArchR)

# 1. Set reference genome build --> reference used for cellranger is hg38
addArchRThreads(20)
addArchRGenome("hg38") # synonym for GRCh38

# 2. Function to read valid barcodes from 10X CellRanger
#' Get Valid Barcodes from 10x Cell Ranger output to pre-filter barcodes
#' 
#' This function will read in processed 10x cell ranger files and identify barcodes that are associated with a cell that passed QC.
#' 
#' @param csvFiles A character vector of names from 10x CSV files to be read in for identification of valid cell barcodes.
#' @param sampleNames A character vector containing the sample names to be associated with each individual entry in `csvFiles`.
#' @export
getValidBarcodes10X <- function(
  csvFiles = NULL, 
  sampleNames = NULL
  ){

  if(length(sampleNames) != length(csvFiles)){
    stop("csvFiles and sampleNames must exist!")
  }

  if(!all(file.exists(csvFiles))){
    stop("Not All csvFiles exists!")
  }

  barcodeList <- lapply(seq_along(csvFiles), function(x){
    df <- suppressPackageStartupMessages(
		suppressMessages(suppressWarnings(readr::read_csv(csvFiles[x])))) 
    if("is__cell_barcode" %ni% colnames(df)){
      stop("is__cell_barcode not in colnames of 10x singlecell.csv file! Are you sure inut is correct?")
    }
    as.character(df[which(paste0(df$is__cell_barcode) ==1),]$barcode)
  }) %>% SimpleList
  names(barcodeList) <- sampleNames

  barcodeList
}

# 3. generate list with fragment files
samples <- read.csv("/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/phenotype_clean.csv", header=TRUE)

inputFiles <- list()
validBarcodes <- list()
for (mb in c(1,2)) {
	for (i in which(samples$Main.Batch == mb)){
		if (mb == 1){
			basepath <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20201104_Anna_Froehlich_10X_ATACseq/"
		} else {
			basepath <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20210316_Anna_Froehlich_10X_ATACseq/"
		}
		file <- paste0(basepath,"01_cellranger/SU",samples$SU.Number[i],"_PFC_ATAC/outs/fragments.tsv.gz")
		inputFiles[paste0("SU",samples$SU.Number[i],"_PFC_ATAC")] <- file

		validBarcodeFile <- paste0(basepath,"01_cellranger/SU",samples$SU.Number[i],"_PFC_ATAC/outs/singlecell.csv")
		validBarcodes[paste0("SU",samples$SU.Number[i],"_PFC_ATAC")] <- validBarcodeFile
	}
}
inputFiles <- unlist(inputFiles)
validBarcodes <- unlist(validBarcodes)

# 4. Get valid barcodes according to 10x CellRanger
vb <- getValidBarcodes10X(csvFiles = validBarcodes, sampleNames = names(validBarcodes))

# 3. create arrow files

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  outputNames = paste0("ArrowFiles/",names(inputFiles)),
  validBarcodes = vb,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags =1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)
