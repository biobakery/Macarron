#' Calculate effect size of metabolic features.
#' 
#' Effect size of a feature is the difference in mean log2 transformed abundances in 
#' test and control/reference groups. For the specified meta.data variable (ptype), effect size 
#' is calculated for all test groups against a specificed reference group.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param mod.assn the output of MACARRoN::findModules(); a dataframe containing feature names, 
#' associate annotations and module assignments
#' @param ptype a reference meta.data of choice. 
#' @param ref a reference level/group in ptype
#' @return my.effsiz
#' 
#' @examples 
#' mod.assn <- findModules(mbx, w, annot="HMDB.ID", chem.info=chem.info, 
#'                         min.mms=10, max.mms=30, brk.mms=5)
#' my.effsiz <- calEffSize(mbx, mod.assn, ptype, ref)
#' 
#' @export


calEffSize <- function(se, mod.assn, ptype, ref){
  
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(assay(se))
  fint <- fint[rownames(mod.assn),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  fint <- t(fint)
  
  #Mean abundance of feature in each group of specified meta.data
  grps <- unique(se[[ptype]])
  message("Calculating mean abundance")
  get.mean <- function(g){
    message(g)
    ind <- ind <- se[[ptype]] == g
    m <- sapply(row.names(mod.assn), function(f) mean(fint[ind,f]))
  }
  all.means <- as.data.frame(do.call(cbind, lapply(grps, get.mean)))
  colnames(all.means) <- grps
  
  #Effect size
  message("Calculating effect size")
  get.es <- function(c){
    message(c)
    es <- all.means[,c]- all.means[,ref]
  }
  test.grps <- setdiff(grps, ref) #test groups 
  my.effsiz <- as.data.frame(do.call(cbind, lapply(test.grps, get.es)))
  colnames(my.effsiz) <- test.grps
  rownames(my.effsiz) <- rownames(mod.assn)
  my.effsiz
  }