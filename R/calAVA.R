#' Calculate abundance versus anchor (AVA) of metabolic features.
#' 
#' AVA of a feature is the ratio of its abundance and the most abundant metabolite in the same module 
#' i.e. the "anchor". Anchor is an annotated/known feature if available or just the most abundant metabolic
#' feature. For every feature, mean abundance in each phenotype or condition is calculated and the maximum is considered
#' for AVA calculation. Singletons are assigned an AVA of 1.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp().
#' @param mod.assn the output of MACARRoN::findMacMod().
#' @param ptype metadata of interest. Default: Column 2 of metadata table. 
#' @param anchor.anno anchor identification. Default: column 2 of annotation dataframe.
#' Note: ptype must be consistent across ava, q-value and effect-size calculations.
#' @return mac.ava
#' 
#' @examples 
#' mac.ava <- calAVA(mbx, mod.assn)
#' 
#' @export

calAVA <- function(se, 
                   mod.assn,
                   ptype = NULL,
                   anchor.anno = NULL)
  {
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(assay(se))
  fint <- fint[rownames(mod.assn),]
  fint <- t(fint)
  
  # Setting the metadata
  if(is.null(ptype)){
    ptype <- names(colData(se))[1]
    message(paste0("Metadata chosen for AVA calculation: ",ptype))
  }else{
    ptype = ptype
    message(paste0("Metadata chosen for AVA calculation: ",ptype))
  }
  grps <- unique(se[[ptype]])
  
  # Identification of the anchor in each module
  anchors <- NULL
  modules <- sort(unique(mod.assn$module))
  
  # Function for anchor and the reference phenotype
  find.anchor <- function(m){
    message(paste0("Finding anchor feature for module: ",m))
    if(dim(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])[1] > 0){
      # known features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])
    }else{
      # all features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m),])
    }
    
    # Mean abundances in phenotypes
    means <- NULL
    for (g in grps){
      ind <- se[[ptype]] == g
      if(length(r) > 1){
        gmean <- apply(fint[ind,r],2,function(x) mean(x, na.rm = TRUE))
        d <- data.frame(m, g, max(gmean), stringsAsFactors = FALSE)
        colnames(d) <- c("module","grp","mma")
      }else{
        d <- data.frame(m, g, mean(fint[ind,r], na.rm = TRUE), stringsAsFactors = FALSE)
        colnames(d) <- c("module","grp","mma")
      }
      means <- rbind(means,d)
    }
    tail(means[order(means[3]),],1)
  }
  f.list <- do.call(rbind, lapply(modules, find.anchor))
  mod.assn$feature <- rownames(mod.assn)
  anchors <- as.data.frame(merge(mod.assn, f.list, by = "module"))
  rownames(anchors) <- anchors$feature
  
  # Abundance versus anchor
  # singletons
  sing.ava <- rep(1, length(anchors[which(anchors$module == 0),1]))
  sing.ava <- as.data.frame(sing.ava)
  rownames(sing.ava) <- rownames(anchors[which(anchors$module == 0),])
  colnames(sing.ava) <- "ava"

  
  # module members
  cal.ava <- function(i){
    message(paste0("Calculating AVA for feature: ",i))
    getMean <- function(g){
      ind <- se[[ptype]] == g
      gmean <- mean(fint[ind,i], na.rm = TRUE)
      gmean
    }
    all.means <- sapply(grps, getMean)
    ava.i <- max(all.means)/as.numeric(as.character(anchors[i,"mma"]))
  }
  members <- rownames(anchors[which(anchors$module != 0),])
  memb.ava <- as.data.frame(do.call(rbind, lapply(members, cal.ava)))
  rownames(memb.ava) <- members
  colnames(memb.ava) <- "ava"
  
  # Combining
  ava <- as.data.frame(rbind(sing.ava, memb.ava))
  ava$module <- mod.assn[row.names(ava),"module"]
  ava <- ava[rownames(mod.assn),]
  ava <- ava[,c(2,1)]
  ava$ava <- round(ava$ava, 6)
  
  # Assign anchors
  anno <- as.data.frame(rowData(se))
  if(is.null(anchor.anno)){
    anchor.anno <- colnames(anno)[2]
  }else{
    anchor.anno <- anchor.anno
  }
  
  ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
  ann.mod <- ann.mod[which(ann.mod > 0)]
  ann.mod <- as.data.frame(ann.mod)
  colnames(ann.mod) <- "module"
  
  assignAnchor <- function(m){
      anchor.feature <- rownames(ava[which(ava$module == m & ava$ava == 1),])
      anchor.name <- as.character(anno[anchor.feature, anchor.anno])
    }
  ann.mod$anchor <- as.character(sapply(ann.mod$module, function(m) assignAnchor(m)))
  ann.mod[ann.mod == "character(0)"] <- ""
  ava$anchor <- as.character(sapply(ava$module, function(m) ann.mod[which(ann.mod$module == m),2]))
  ava[ava == "character(0)"] <- ""
  ava
}