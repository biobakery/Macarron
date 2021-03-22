#' Calculate relative abundance of metabolic features.
#' 
#' Relative abundance of a feature is the ratio of its abundance and the most abundant metabolite in the same module 
#' i.e. the "anchor". Anchor is an annotated/known feature if available or just the most abundant metabolic
#' feature. For every feature, mean abundance in each condition is calculated and the maximum is considered
#' for relative abundance calculation. Singletons are assigned a relative abundance of 1.
#' 
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param mod.assn the output of MACARRoN::findModules(); a dataframe containing feature names, 
#' associate annotations and module assignments
#' @param ptype a reference meta.data of choice. 
#' @return my.relab 
#' 
#' @examples 
#' mod.assn <- findModules(mbx, w, annot="HMDB.ID", chem.info=chem.info, 
#'                         min.mms=10, max.mms=30, brk.mms=5)
#' my.relab <- calRelAbn(mbx, mod.assn, ptype)
#' 
#' @export

calRelAbn <- function(se, mod.assn, ptype){
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(assay(se))
  fint <- fint[rownames(mod.assn),]
  fint <- t(fint)
  
  # Identification of the anchor in each module
  anchors <- NULL
  modules <- sort(unique(mod.assn$module))
  grps <- unique(se[[ptype]])
  
  find.anchor <- function(m){
    message(paste0("Finding anchor feature for module ",m))
    if(dim(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])[1] > 0){
      #known features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])
    }else{
      #all features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m),])
    }
    
    means <- NULL
    for (g in grps){
      ind <- se[[ptype]] == g
      if(length(r) > 1){
        gmean <- apply(fint[ind,r],2,function(x) mean(x, na.rm=TRUE))
        d <- data.frame(m, g, max(gmean), stringsAsFactors = FALSE)
        colnames(d) <- c("module","grp","mma")
      }else{
        d <- data.frame(m, g, mean(fint[ind,r], na.rm=TRUE), stringsAsFactors = FALSE)
        colnames(d) <- c("module","grp","mma")
      }
      means <- rbind(means,d)
    }
    tail(means[order(means[3]),],1)
  }
  f.list <- do.call(rbind, lapply(modules, find.anchor))
  mod.assn$feature <- rownames(mod.assn)
  anchors <- merge(mod.assn, f.list, by="module")
  rownames(anchors) <- anchors$feature
  anchors
  
  #Relative abundance
  #singletons
  sing.ra <- rep(1, length(anchors[which(anchors$module==0),1]))
  sing.ra <- as.data.frame(sing.ra)
  rownames(sing.ra) <- rownames(anchors[which(anchors$module==0),])
  colnames(sing.ra) <- "rel.abun."
  
  #module members
  cal.relab <- function(i){
    message(paste0("Calculating relative abundance for feature ",i))
    getMean <- function(g){
      ind <- se[[ptype]] == g
      gmean <- mean(fint[ind,i], na.rm=TRUE)
      gmean
    }
    all.means <- sapply(grps, getMean)
    relab.i <- max(all.means)/as.numeric(as.character(anchors[i,"mma"]))
  }
  members <- rownames(anchors[which(anchors$module!=0),])
  memb.ra <- as.data.frame(do.call(rbind, lapply(members, cal.relab)))
  rownames(memb.ra) <- members
  colnames(memb.ra) <- "rel.abun."
  
  #Combining 
  my.relab <- as.data.frame(rbind(sing.ra, memb.ra))
  my.relab$module <- mod.assn[row.names(my.relab),"module"]
  my.relab <- my.relab[rownames(mod.assn),]
  my.relab
}
