## functions to detect and remove outliers and to remove batch effects

# Outlier detection on pseudosamples of one cell type
outlier_removal <- function(dds, prefix) {
  outlier_bool <- TRUE
  iter <- 1
  while (outlier_bool) {
    # Run PCA
    dds <- estimateSizeFactors(dds)
    vsd <- vst(dds, blind = TRUE)
    pca_result <- pca(vsd@assays@data[[1]])
    pc <-
      as.data.frame(pca_result$rotated)
    pc$sample <- rownames(pc)
    
    # Label any sample which is more than 3SD away from the mean in PC1 as outlier
    outlier <-
      which(abs(pc[, 1] - mean(pc[, 1])) > (3 * sd(pc[, 1])))
    pc$outlier <- FALSE
    pc$outlier[outlier] <- TRUE
    
    ggplot(pc, aes(
      x = PC1,
      y = PC2,
      color = outlier,
      label = sample
    )) +
      geom_point(size = 3) +
      geom_text_repel() +
      xlab(paste0("PC1: ", pca_result$variance[1])) +
      ylab(paste0("PC2: ", pca_result$variance[2])) +
      coord_fixed() +
      ggtitle("Principal Component Analysis")
    ggsave(paste0(prefix, "_outlier_iter", iter, ".png"))
    
    # Subset deseq object / remove outliers
    keep_names <-
      setdiff(colData(dds)$sample_id, pc$sample[outlier])
    dds <- dds[, keep_names]
    iter <- iter + 1
    if (length(outlier) == 0) {
      outlier_bool <- FALSE
    }
  }
  
  # Rerun normalization etc on outlier filtered dataset
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)
  pca_result <- pca(vsd@assays@data[[1]])
  pc <-
    as.data.frame(pca_result$rotated)
  pc$Status <- colData(dds)$Status
  
  ggplot(pc, aes(x = PC1, y = PC2, color = Status)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", pca_result$variance[1])) +
    ylab(paste0("PC2: ", pca_result$variance[2])) +
    coord_fixed() +
    ggtitle("PCA: Outlier Deleted")
  ggsave(paste0(prefix, "_outlierDeleted.png"))
  
  png(
    filename = paste0(prefix, "_pairsplot_outlierDeleted.png"),
    width = 900,
    height = 900
  )
  print(splom(
    as.data.frame(pc[, 1:10]),
    col = colData(dds)$Status,
    cex = 2,
    pch = '*'
  ))
  dev.off()
  
  # # add new PCs to colData
  # colData(dds)[, c((ncol(colData(dds))+1):(ncol(colData(dds))+9))] <- pc[, c(1:10)]
  
  return(dds)
}


# Batch identification with SVA
batch_sva <- function(dds, mod_dds){
  
  # 1. Surrogate Variable Analysis (SVA) ----
  mod <- model.matrix(as.formula(mod_dds), colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  # Calculate SVs (on normalized counts)
  # Comment: don't use num.sv function, apparently it's only for microarray data
  norm <- counts(dds, normalized = TRUE)
  svobj <- svaseq(norm, mod, mod0)
  n.sv <- svobj$n.sv #Number of significant surrogate variables
  # Add significant SVs to covariates
  coln <- colnames(colData(dds))
  cov_pc <- cbind(colData(dds),svobj$sv[,1:n.sv])
  colnames(cov_pc) <- c(coln, paste0("SV", seq(1:n.sv)))
  
  
  # 2. Canonical Correlation Analysis (CCA) ----
  form <- as.formula(paste0(des,"+",
                    paste(paste0("SV", seq(1:n.sv)), collapse = "+")))
  
  # Calculate the correlation coefficients
  C <- canCorPairs(form, cov_pc)
  # Plot the results using Canonical correlation
  png(filename = paste0(basedir, "plots/cca_sorted", cluster, "_", des, ".png"), width = 800, height = 800)
  plotCorrMatrix(C)
  dev.off()
  png(filename = paste0(basedir, "plots/cca_unsorted", cluster, "_", des, ".png"), width = 800, height = 800)
  plotCorrMatrix(C, sort = FALSE)
  dev.off()
  
  return(list("cov_pc" = cov_pc, "n.sv" = n.sv))
}


# Define own method that is able to remove 3 categorical batch variables

remove3BatchEffects <- function(x,batch=NULL,batch2=NULL,batch3=NULL,
                                covariates=NULL,design=matrix(1,ncol(x),1),...)
  #  Remove batch effects from matrix of expression data
  #  Adapted from Gordon Smyth and Carolyn de Graaf
  #  Created 1 Aug 2008. Last revised 1 June 2014.
{
  if(is.null(batch) && is.null(batch2) && is.null(batch3) && is.null(covariates)) return(as.matrix(x))
  if(!is.null(batch)) {
    batch <- as.factor(batch)
    contrasts(batch) <- contr.sum(levels(batch))
    batch <- model.matrix(~batch)[,-1,drop=FALSE]
  }
  if(!is.null(batch2)) {
    batch2 <- as.factor(batch2)
    contrasts(batch2) <- contr.sum(levels(batch2))
    batch2 <- model.matrix(~batch2)[,-1,drop=FALSE]
  }
  if(!is.null(batch3)) {
    batch3 <- as.factor(batch3)
    contrasts(batch3) <- contr.sum(levels(batch3))
    batch3 <- model.matrix(~batch3)[,-1,drop=FALSE]
  }
  if(!is.null(covariates)) covariates <- as.matrix(covariates)
  X.batch <- cbind(batch,batch2,batch3,covariates)
  fit <- lmFit(x,cbind(design,X.batch),...)
  beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
  beta[is.na(beta)] <- 0
  as.matrix(x) - beta %*% t(X.batch)
}
