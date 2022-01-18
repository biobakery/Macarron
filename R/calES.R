#' Calculate effect size of differential abundance of metabolic features.
#' 
#' Effect size of a metabolic feature is the difference in mean log2 transformed abundances in 
#' test and control (reference) samples. For the specified metadata variable, effect size 
#' is calculated for all test categories against the reference category.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp().
#' @param mac.qval the output of MACARRoN::calQval().
#' @return mac.efs
#' 
#' @examples 
#' mac.es <- calES(se, mac.qval)
#' 
#' @export

calES <- function(se,
                  mac.qval)
{
  # Abundance matrix
  fint <- as.data.frame(SummarizedExperiment::assay(se))
  fint <- fint[unique(mac.qval$feature),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  fint <- t(fint)
  
  # Get phenotype from mac.qval
  ptype <- unique(mac.qval$metadata)
  grps <- unique(se[[ptype]])
  
  # Mean abundance of feature in each group of ptype
  get.mean <- function(g){
    message(paste0("Calculating mean abundance in: ",g))
    ind <- se[[ptype]] == g
    m <- sapply(colnames(fint), function(f) mean(fint[ind,f]))
  }
  all.means <- as.data.frame(do.call(cbind, lapply(grps, get.mean)))
  colnames(all.means) <- grps
  
  # Effect size
  test.grps <- unique(mac.qval$value)
  ref.grp <- setdiff(grps, test.grps)
  get.es <- function(tg){
    message(paste0("Calculating effect size in: ",tg))
    es <- all.means[,tg]- all.means[,ref.grp]
  }
  mac.es <- as.data.frame(do.call(cbind, lapply(test.grps, get.es)))
  colnames(mac.es) <- test.grps
  rownames(mac.es) <- colnames(fint)
  mac.es
}