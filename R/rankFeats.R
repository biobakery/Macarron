#' Rank metabolic features and prioritize based on predicted bioactivity
#' 
#' Metabolic features are ranked based on relative abundance, effect size and q-value of
#' differential abundance. The harmonic mean of these three ranks is calculated and used as 
#' the meta-rank to prioritize potentially bioactive features in a user-specified 
#' environment or phenotype (a level in a meta.data variable). Top-ranked features have
#' good relative abundance, and are significantly perturbed in the specified environment/
#' phenotype.
#' 
#' @param my.relab relative abundance of features calculated using MACARRoN::calRelAbun()
#' @param my.effsiz effect size of differential abundance of features in the specified phenotype/
#' environment calculated using MACARRoN::calEffSize()
#' @param my.qvalue q-value of differential abundance of features in the specified phenotype/
#' environment calculated using MACARRoN::calQvalue()
#' @param mod.assn the output of MACARRoN::findModules(); a dataframe containing feature names, 
#' associate annotations and module assignments
#' @param category the specified environment/phenotype
#' @return ranked.features
#' 
#' @examples
#' mod.assn <- findModules(mbx, w, annot="HMDB.ID", chem.info=chem.info, 
#'                         min.mms=10, max.mms=30, brk.mms=5)
#' my.relab <- calRelAbn(mbx, mod.assn, ptype)
#' my.effsiz <- calEffSize(mbx, mod.assn, ptype, ref)
#' my.qvalue <- calQvalue(mbx, mod.assn, ptype, ref, fix.eff, ran.eff, folder.name, ncores)
#' ranked.features <- rankFeats(my.relab, my.effsiz, my.qvalue, mod.assn, category=category)
#' 
#' @export


rankFeats <- function(se, my.relab, my.effsiz, my.qvalue, mod.assn, category, annot){
  # packages
  library(psych)
  
  # collecting all parameters
  message("collecting all features")
  qval <- my.qvalue[which(my.qvalue$value == category),]
  rownames(qval) <- qval$feature
  qval <- qval[rownames(my.effsiz),]
  all.param <- as.data.frame(cbind(my.relab[,"rel.abun."], my.effsiz[,category], qval[,"qvalue"]))
  colnames(all.param) <- c("relab","efs","qval")
  rownames(all.param) <- rownames(my.effsiz)
  
  # direction of perturbation and annotations
  message("assigning direction of perturbation")
  all.param$status <- ""
  all.param[which(all.param$efs < 0),"status"] <- paste0("depleted in ", category)
  all.param[which(all.param$efs > 0),"status"] <- paste0("enriched in ", category)
  all.param$efs <- abs(all.param$efs)
  all.param$module <- mod.assn[row.names(all.param),2]
  all.param$annotation <- mod.assn[row.names(all.param),1]
  all.param$association <- 0
  all.param[which(all.param$annotation != ""),"association"] <- 1
  
  mods.with.stds <- unique(all.param[which(all.param$annotation != "" & all.param$module != 0),"module"])
  all.param[which(all.param$module %in% mods.with.stds),"association"] <- 1
  
  
  # ranks
  message("assigning ranks")
  all.param$relab_rank <- rank(all.param$relab)
  all.param$efs_rank <- rank(all.param$efs)
  all.param$qval_rank <- rank(-all.param$qval)
  all.param$assoc_rank <- rank(all.param$association)
  
  # meta rank
  message("calculating meta-rank and prioritizing features")
  all.param$meta_rank <- harmonic.mean(t(all.param[,8:11]))
  ranked.features <- all.param[order(-all.param$meta_rank),]
  ranked.features$priority <- seq(1:nrow(ranked.features))
  ranked.features <- as.data.frame(cbind(ranked.features, rowData(se)[rownames(ranked.features),annot]))
  ranked.features
}