## Functions to enable TF motif enrichment analysis

.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}


.computeEnrichment <- function(matches = NULL, compare = NULL, background = NULL){

  matches <- .getAssay(matches,  grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  #Enrichment
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
     pOut$BackgroundFrequency[x], # Number of all successes in background
     pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
     pOut$nCompare[x], # Number that were drawn
     lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)

  #Minus Log10 Padj
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

  pOut

}

#' Peak Annotation Hypergeometric Enrichment in a given set of Peaks.
#' 
#' This function will perform hypergeometric enrichment of a given peak annotation within the given peaks of interest.
#' 
#' @param seMarker  A `GRanges` object with peaks of interest. Name kept for ease of adaptations.
#' @param ArchRProj An `ArchRProject` object.
#' @param peakAnnotation A `peakAnnotation` object in the provided `ArchRProject` to be used for hypergeometric testing.
#' @param matches A custom `peakAnnotation` matches object used as input for the hypergeometric test. See
#' `motifmatchr::matchmotifs()` for additional information.
#' @param background A string that indicates whether to use a background set of matched peaks to compare against ("bgdPeaks") or all peaks ("all").
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
peakAnnoEnrichmentAdapted <- function(
  seMarker = NULL,
  ArchRProj = NULL,
  peakAnnotation = NULL,
  matches = NULL,
  background = "all",
  logFile = createLogFile("peakAnnoEnrichment")
  ){

#   .validInput(input = seMarker, name = "seMarker", valid = c("GRanges"))
#   .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
#   .validInput(input = peakAnnotation, name = "peakAnnotation", valid = c("character", "null"))
#   .validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", "null"))
#   .validInput(input = background, name = "background", valid = c("character"))
#   .validInput(input = logFile, name = "logFile", valid = c("character"))

#   tstart <- Sys.time()
#   .startLogging(logFile = logFile)
#   .logThis(mget(names(formals()),sys.frame(sys.nframe())), "peakAnnoEnrichment Input-Parameters", logFile = logFile)
#
#   if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
#     stop("Only markers identified from PeakMatrix can be used!")
#   }

  if(is.null(matches)){
    matches <- getMatches(ArchRProj, peakAnnotation)
  }
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  mcols(r1) <- NULL

  r2 <- getPeakSet(ArchRProj)
  pr2 <- paste(seqnames(r2),start(r2),end(r2),sep="_")
  mcols(r2) <- NULL

# r3 <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, rowData(seMarker)$end))
  r3 <- seMarker
  pr3 <- paste(seqnames(r3),start(r3),end(r3),sep="_")
  mcols(r3) <- NULL

#   .logThis(r1, "Peaks-Matches", logFile = logFile)
#   .logThis(r2, "Peaks-ArchRProj", logFile = logFile)
#   .logThis(r3, "Peaks-SeMarker", logFile = logFile)

#   .logThis(pr1, "Peaks-Pasted-Matches", logFile = logFile)
#   .logThis(pr2, "Peaks-Pasted-ArchRProj", logFile = logFile)
#   .logThis(pr3, "Peaks-Pasted-SeMarker", logFile = logFile)

  if(length(which(pr1 %ni% pr2)) != 0){
    stop("Peaks from matches do not match peakSet in ArchRProj!")
  }

#   if(length(which(pr2 %ni% pr3)) != 0){
#     stop("Peaks from seMarker do not match peakSet in ArchRProj!")
#   }

  rownames(matches) <- pr1
  #matches <- matches[pr3, ]

  #Evaluate AssayNames
#   assayNames <- names(SummarizedExperiment::assays(seMarker))
#   for(an in assayNames){
#     eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
#   }
#   passMat <- eval(parse(text=cutOff))
#   passMat[is.na(passMat)] <- FALSE
#   for(an in assayNames){
#     eval(parse(text=paste0("rm(",an,")")))
#   }

  # TODO: fix this, as currently only the all option really works
  if(tolower(background) %in% c("backgroundpeaks", "bgdpeaks", "background", "bgd")){
    method <- "bgd"
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  }else{
    method <- "all"
  }

#   enrichList <- lapply(seq_len(ncol(seMarker)), function(x){
#     #.logDiffTime(sprintf("Computing Enrichments %s of %s", x, ncol(seMarker)), t1 = tstart, verbose = TRUE, logFile = logFile)
#     #idx <- which(passMat[, x])
#     if(method == "bgd"){
#       .computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx,])))
#     }else{
#       .computeEnrichment(matches, idx, seq_len(nrow(matches)))
#     }
#   }) %>% SimpleList
#   names(enrichList) <- colnames(seMarker)

  #idx <- which(rowRanges(matches)$idx %in% seMarker$idx)
  idx <- which(rownames(matches) %in% pr3)
  enrichment <- .computeEnrichment(matches, idx, seq_len(nrow(matches)))

#   assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
#     d <- lapply(seq_along(enrichList), function(y){
#       enrichList[[y]][colnames(matches),x,drop=FALSE]
#     }) %>% Reduce("cbind",.)
#     colnames(d) <- names(enrichList)
#     d
#   }) %>% SimpleList
#   names(assays) <- colnames(enrichList[[1]])
#   assays <- rev(assays)
#   out <- SummarizedExperiment::SummarizedExperiment(assays=assays)

#   .endLogging(logFile = logFile)

#   out

enrichment

}


fdr_perCelltype <- function(data_df, pvalue_column){

	data_list <- split(data_df, f = data_df$celltype)
	# add column with corrected p-values for each celltype
	data_list <- lapply(data_list, function(x) x %>% dplyr::mutate(FDR = p.adjust(x[[pvalue_column]], method="fdr")))
	out_df <- bind_rows(data_list)
	#out_df$celltype <- factor(out_df$celltype, levels=celltypes)
	out_df$sig <- ifelse(out_df$FDR <= 0.1, TRUE, FALSE)

	return(out_df)
}