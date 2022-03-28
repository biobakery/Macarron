#' Create a biweight midcorrelation (WGCNA::bicor()) based distance matrix.
#' 
#' @param se SummarizedExperiment object created using Macarron::prepInput().
#' @param metadata_variable metadata column identifying phenotypes/conditions to be used to evaluate prevalence of features. Default = Column 1 of metadata dataframe.
#' @param min_prevalence prevalence threshold (percentage). Default = 0.7.
#' @param execution_mode "serial" or "multi" processing with BiocParallel. Default: "serial" (recommended for laptops). 
#' "multi" may be used when running Macarron on a cluster. 
#' @param optimize.for runtime or memory.
#' 
#' Features present (i.e. not NA) in "min_prevalence" of samples in each category of a "metadata_variable" will be considered 
#' e.g. if min_prevalence is 0.7 and metadata_variable has 2 categories A and B, union of (i) features present in at least 70% of A samples
#' and (ii) features present in at least 70% of B samples, will be considered for distance matrix generation. 
#' Correlation between feature abundances are is calculated using WGCNA::bicor().
#' @return w distance matrix where distance = 1-bicor^3
#' 
#' @examples 
#' prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
#' abundances_df = read.csv(file = prism_abundances, row.names = 1)
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
#' metadata_df = read.csv(file = prism_metadata, row.names = 1)
#' mbx <- Macarron::prepInput(input_abundances = abundances_df,
#'                             input_annotations = annotations_df,
#'                             input_metadata = metadata_df)
#' w <- Macarron::makeDisMat(se = mbx)
#' 
#' 
#' @export

makeDisMat <- function(se, 
                       metadata_variable = 1,
                       min_prevalence = 0.7,
                       execution_mode = "serial",
                       optimize.for = c("runtime", "memory"))
{
  optimize.for <- match.arg(optimize.for)
  opt.mem <- optimize.for == "memory"
  
  # Abundance matrix
  mat <- DelayedArray::DelayedArray(SummarizedExperiment::assay(se))
  
  # Phenotype i.e. groups/conditions
  if(is.character(metadata_variable)){
    metadata_variable <- metadata_variable
  }else{
    metadata_variable <- names(SummarizedExperiment::colData(se))[metadata_variable]
  }
  phenotypes <- unique(se[[metadata_variable]])
  
  # Function: Get features that satisfy prevalence threshold in each condition
  .getIds <- function(phenotype)
  {
    ind <- se[[metadata_variable]] == phenotype
    smat <- mat[,ind]
    ind <- DelayedArray::rowMeans(is.na(smat)) <= 1 - min_prevalence
    ind
  }
  
  # Apply on all groups/conditions
  ind <- vapply(phenotypes, .getIds, logical(nrow(mat)))
  
  # Union of features
  ind <- apply(ind,1,any)
  mat <- mat[ind,]
  message(paste0(nrow(mat), " features pass chosen minimum prevalence threshold of ", min_prevalence, "."))
  
  # Compute correlation and distance matrices
  if(nrow(mat) <= 25000){
    .getCorMat <- function(phenotype)
    {
      message(paste0("Calculating pairwise correlations in phenotype: ", phenotype))
      inx <- se[[metadata_variable]] == phenotype
      ind <- DelayedArray::rowMeans(!is.na(mat[,inx]))>=min_prevalence
      gmat <- mat[ind,inx]
      tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
      rownames(tmp) <- rownames(mat)
      colnames(tmp) <- colnames(gmat)
      tmp <- DelayedArray::DelayedArray(tmp)
      tmp[rownames(gmat), colnames(gmat)] <- gmat
      tmp <- log2(tmp)
      options(warn=-1)
      cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs", quick=0.05)
      options(warn=0)
      cmat[is.na(cmat)] <- 0
      cmat[cmat < 0] <- 0
      if(opt.mem) cmat <- as(cmat, "dsyMatrix")
      cmat
    }
    if(execution_mode == "serial"){
      exe.choice <- BiocParallel::SerialParam()
    }else if(execution_mode == "multi"){
      exe.choice <- BiocParallel::MulticoreParam()
    }
    # Apply on all groups/conditions
    cmats <- BiocParallel::bplapply(phenotypes, .getCorMat, BPPARAM = exe.choice)
    # Keep the best observed positive correlation for each pair of features
    mmat <- do.call(pmax, c(cmats, na.rm=TRUE))
    }else{
      for (phenotype in phenotypes){
        message(paste0("Calculating pairwise correlations in phenotype: ", phenotype))
        inx <- se[[metadata_variable]] == phenotype
        ind <- DelayedArray::rowMeans(!is.na(mat[,inx]))>=min_prevalence
        gmat <- mat[ind,inx]
        tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
        rownames(tmp) <- rownames(mat)
        colnames(tmp) <- colnames(gmat)
        tmp <- DelayedArray::DelayedArray(tmp)
        tmp[rownames(gmat), colnames(gmat)] <- gmat
        tmp <- log2(tmp)
        options(warn=-1)
        cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs", quick=0.05)
        options(warn=0)
        cmat[is.na(cmat)] <- 0
        cmat[cmat < 0] <- 0
        assign(paste0(phenotype,"_ff"),ff::as.ff(cmat))
        rm(cmat)
       }
    # Keep the best observed positive correlation for each pair of features
    mmat <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
    for(j in ls(pattern="_ff")){
      for(k in 1:nrow(mmat)){
        mmat[,k] <- pmax(mmat[,k],get(j)[,k],na.rm=TRUE)
        }
    }
    rownames(mmat) <- rownames(mat)
    colnames(mmat) <- rownames(mat)
    }
  
  # Beta-scaling
  mmat = mmat^3 
    
  # Distance matrix
  w = 1 - mmat
  message(paste0("Distance matrix with ",nrow(w)," features created."))
  w
}
  
