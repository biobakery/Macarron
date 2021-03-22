#' Create a biweight midcorrelation (WGCNA::bicor()) based distance matrix
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param ptype meta.data (phenotype/condition) to be used to evaluate prevalence of features.
#' @param preval prevalence threshold (percentage). Default is 0.7
#' @param BPPARAM serial or parallel processing with BiocParallel. Default: SerialParam (recommended for laptops). MulticoreParam() can be used when running MACARRoN on a cluster. 
#' Features present (i.e. not NA) in "preval" of samples in each category of a "ptype" will be considered 
#' e.g. if preval is 0.7 and ptype has 2 categories A and B, union of (i) features present in 70% of A samples
#' and (ii) features present in 70% of B samples will be considered for distance matrix generation. 
#' Correlation between feature abundances are is calculated using WGCNA::bicor()
#' @return w distance matrix where distance is 1-bicor^3
#' 
#' @examples 
#' mbx <- makeSumExp(feat_int, feat_anno, exp_meta)
#' w <- makeDisMat(mbx, ptype="diagnosis", preval=0.7)
#' 
#' @export


makeDisMat <- function(se, ptype, preval=0.7,
                         BPPARAM = BiocParallel::SerialParam(),
                         optimize.for = c("runtime", "memory"))
{
  optimize.for <- match.arg(optimize.for)
  opt.mem <- optimize.for == "memory"
  
  # packages
  library(BiocParallel)
  library(DelayedArray)
  library(WGCNA)
  library(ff)
  library(ffbase)
  library(data.table)
  
  # Abundance matrix
  mat <- DelayedArray::DelayedArray(SummarizedExperiment::assay(se))
  
  # Phenotype i.e. groups/conditions
  grps <- unique(se[[ptype]])
  
  # Function: Get features that satisfy prevalence threshold in each condition
  .getIds <- function(g)
  {
    ind <- se[[ptype]] == g
    smat <- mat[,ind]
    ind <- rowMeans(is.na(smat)) <= 1 - preval
    ind
  }
  
  # Apply on all groups/conditions
  ind <- vapply(grps, .getIds, logical(nrow(mat)))
  
  # Union of features
  ind <- apply(ind,1,any)
  mat <- mat[ind,]
  message(paste0(nrow(mat), " features are selected."))
  
  
  # Compute correlation and distance matrices
  # Small datasets
  if(nrow(mat) <= 30000){
    
    #Function: Correlation matrix
    .getCorMat <- function(g)
    {
      message(g)
      inx <- se[[ptype]] == g
      gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=preval),inx]
      tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
      rownames(tmp) <- rownames(mat)
      colnames(tmp) <- colnames(gmat)
      tmp <- DelayedArray::DelayedArray(tmp)
      tmp[rownames(gmat), colnames(gmat)] <- gmat
      tmp[is.na(tmp)] <- 0
      tmp <- log2(tmp + 1)
      cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs")
      cmat[is.na(cmat)] <- 0
      cmat[cmat < 0] <- 0
      if(opt.mem) cmat <- as(cmat, "dsyMatrix")
      cmat
    }
    
    # Apply on all groups/conditions
    cmats <- BiocParallel::bplapply(grps, .getCorMat, BPPARAM = BPPARAM)
    
    # Keep the best observed positive correlation for each pair of features
    mmat <- do.call(pmax, c(cmats, na.rm=TRUE))
    
    # Beta-scaling to ensure power law distribution
    mmat <- mmat^3
    
    # Distance matrix
    w = 1 - mmat
    message(paste0("Distance matrix with ",nrow(w)," features created."))
    w
    }else
      # Large datasets
      {for (g in grps){
        message(g)
        inx <- se[[ptype]] == g
        gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=preval),inx]
        tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
        rownames(tmp) <- rownames(mat)
        colnames(tmp) <- colnames(gmat)
        tmp <- DelayedArray::DelayedArray(tmp)
        tmp[rownames(gmat), colnames(gmat)] <- gmat
        tmp[is.na(tmp)] <- 0
        tmp <- log2(tmp + 1)
        cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs")
        message(paste0(g," cmat created"))
        assign(paste0(g,"_ff"),as.ff(cmat))
        message(paste0(g,"_ff created"))
        rm(cmat)
        }
    
    # Keep the best observed positive correlation for each pair of features
    mmat <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
    for(j in ls(pattern="_ff")){
      for(k in 1:nrow(mmat)){
        mmat[,k] <- pmax(mmat[,k],get(j)[,k],na.rm=TRUE)
        }
      message(paste0(j," compared"))
    }
    rownames(mmat) <- rownames(mat)
    colnames(mmat) <- rownames(mat)
    mmat = mmat^3 
    
    # Distance matrix
    w = 1 - mmat
    message(paste0("Distance matrix with ",nrow(w)," features created."))
    w}
}
  
