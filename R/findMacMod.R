#' Cluster metabolic features based on covarying abundances into modules
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp().
#' @param w distance matrix from function MACARRoN::makeDisMat().
#' @param chem_tax chemical taxonomy file with 3 columns specifying annotation, subclass and class of annotated features. 
#' Can be created using the decorateID.R utility of MACARRoN. 
#' Annotation specified with "annot" and annotation in the first column of the chemical taxonomy file must match.
#' @param annot a feature annotation of choice e.g. HMDB ID/Pubchem CID/Metabolite name. Default: Column 1 in annotation table.
#' @param mms minimum module size to be used for module identification with dynamicTreeCut::cutreeDynamic(). 
#' Default is cube root of number of prevalent features.
#' @param evaluateMOS examine measure of success for modules identified using mms, mms + 5, mms + 10, mms - 5, mms - 10
#' 
#' @return mod.assn metabolic features clustered into "modules" based on covarying abundances
#' 
#' @examples 
#' mod.assn <- findMacMod(se, w, chem_tax)
#' 
#' @importFrom dynamicTreeCut cutreeDynamic
#' 
#' @export

findMacMod <- function(se, 
                        w, 
                        chem_tax = NULL,
                        annot = NULL,
                        mms = NULL,
                        evaluateMOS = TRUE)
{
  # packages
  requireNamespace("dynamicTreeCut", quietly = TRUE)
  
  # Construct tree 
  tree <- hclust(as.dist(w), method="average")
  message("Tree constructed")
  
  # Setting the minimum module size
  if(is.null(mms)){
    mms = round((nrow(w))^(1/3))
  }else{
    mms = as.numeric(as.character(mms))
  }
  
  # Module assignments
  anno <- as.data.frame(SummarizedExperiment::rowData(se))
  if(is.null(annot)){
    mod.assn <- as.data.frame(anno[colnames(w),1])
  }else{
    mod.assn <- as.data.frame(anno[colnames(w),annot])
  }
  rownames(mod.assn) <- colnames(w)
  
  #-------------------------------------
  # Measures of success
  #-------------------------------------
  if(isTRUE(evaluateMOS)){
    
    message("Evaluating measures of success")
    
    # Range of mms for measures of success
    mms.list <- c(mms - 10, mms - 5, mms, mms + 5, mms + 10)
    mms.list <- mms.list[which(mms.list > 1)]
    
    # Measures of success
    sing <- NULL # singletons
    totc <- NULL # total modules
    pann <- NULL # % annotated modules
    hscc <- NULL # % features with the same annotation that are assigned to the same module
    maxc <- NULL # max classes per module
    maxs <- NULL # max subclasses per module
    perc <- NULL # 90th percentile classes per module
    pers <- NULL # 90th percentile subclasses per module
    fham <- NULL # % features in homogeneously annotated modules
    
    for (n in mms.list){
      mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                              distM = as.matrix(w), 
                                                              deepSplit = TRUE, 
                                                              pamRespectsDendro = TRUE,
                                                              minClusterSize = n))
      # singletons
      sing <- rbind(sing, length(which(mod.assn[,2] == 0)))
      
      # total modules
      totc <- rbind(totc, max(mod.assn[,2]))
      
      # % annotated modules
      ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
      ann.mod <- ann.mod[which(ann.mod > 0)]
      pann <- rbind(pann, round((length(ann.mod)*100)/max(mod.assn[,2]),2))
      
      # % features with the same annotation that are assigned to the same module
      same.mod <- NULL
      diff.mod <- NULL
      times.obs <- as.data.frame(table(mod.assn[,1][which(mod.assn[,1] != "")]))
      times.obs <- times.obs[which(times.obs$Freq > 1),]
      for (t in times.obs$Var1){
        if(length(unique(mod.assn[which(mod.assn[,1] == t),2])) > 1){
          diff.mod <- rbind(diff.mod,t)}else{
            same.mod <- rbind(same.mod,t)
          }
      }
      hscc <- rbind(hscc,round(length(same.mod)*100/(length(same.mod)+length(diff.mod)),2))
      
      # chemical homogeneity of modules
      cls <- NULL
      scl <- NULL
      mod.assn.ann <- mod.assn[which(mod.assn[,2] %in% ann.mod & mod.assn[,1] != ""),]
      for (i in unique(mod.assn.ann[,2])){
        dat <- mod.assn.ann[which(mod.assn.ann[,2] == i & mod.assn.ann[,1] != ""),]
        mod.tax <- NULL
        for (d in unique(dat[,1])){
          feat.tax <- chem_tax[which(chem_tax[,1] == d),]
          mod.tax <- rbind(mod.tax,feat.tax)
        }
        cls <- rbind(scl, length(unique(mod.tax[,3])))
        scl <- rbind(scl, length(unique(mod.tax[,2])))
      }
      maxs <- rbind(maxs, max(scl))
      maxc <- rbind(maxc, max(cls))
      pers <- rbind(pers, quantile(scl,0.9))
      perc <- rbind(perc, quantile(cls,0.9))
      
      # % features in homogeneously annotated modules
      dat <- as.data.frame(unique(mod.assn.ann))
      dat$class <- sapply(as.character(dat[,1]), 
                          function(x) as.character(chem_tax[which(chem_tax[,1] == x),3]))
      rownames(dat) <- NULL
      dat$class <- as.character(dat$class)
      dat$class[which(dat$class == "character(0)")] <- ""
      dat <- dat[which(dat$class != ""),]
      mods.with.tax <- as.data.frame(sort(unique(dat[,2])))
      names(mods.with.tax) <- "module"
      
      # function to estimate homogeneity of each module
      findClassHomo <- function(m){
        df <- dat[which(dat[,2] == m),]
        annotated.features <- nrow(df)
        class.counts <- as.data.frame(table(df$class))
        max(class.counts$Freq)/annotated.features
      }
      
      mods.with.tax$homo <- sapply(mods.with.tax$module, function(m) findClassHomo(m))
      mods.with.tax$homo <- as.numeric(as.character(mods.with.tax$homo))
      homogeneous.mods <- unique(mods.with.tax[which(mods.with.tax$homo >= 0.75),"module"])
      
      perc.feats <- round((nrow(mod.assn[which(mod.assn[,2] %in% homogeneous.mods),])*100)/nrow(mod.assn),2)
      fham <- rbind(fham, perc.feats)
    }
    
    # Writing results to file
    mac.mos <- data.frame(cbind(mms.list, totc, sing, pann, hscc, maxc, perc, maxs, pers, fham))
    rownames(mac.mos) <- NULL
    colnames(mac.mos) <- c("Minimum module size (MMS)",
                           "Total modules",
                           "Singletons",
                           "% Annotated modules",
                           "% Consistent assignments",
                           "Max classes per module",
                           "90p classes per module",
                           "Max subclasses per module",
                           "90p subclasses per module",
                           "% Features in HAM")
                      
    write.csv(mac.mos, file="MAC_modules_measures_of_success.csv")
  }
  
  #-------------------------------------
  # Module assignments
  #-------------------------------------  
  mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                          distM = as.matrix(w), 
                                                          deepSplit = TRUE, 
                                                          pamRespectsDendro = TRUE,
                                                          minClusterSize = mms))
  if(is.null(annot)){
    colnames(mod.assn) <- c(names(anno)[1],"module")
  }else{
    colnames(mod.assn) <- c(annot,"module")
  }
  
  ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
  ann.mod <- ann.mod[which(ann.mod > 0)]
  ann.mod <- as.data.frame(ann.mod)
  colnames(ann.mod) <- "module"
  
  # Annotate modules with chemical composition information based on annotated features
  assignChemTax <- function(m){
    if(m > 0){
      annotated.features <- unique(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),1])
      classes <- toString(unique(chem_tax[which(chem_tax[,1] %in% annotated.features),3]))
    }else{
      classes <- ""
    }
  }
  ann.mod$classes <- as.character(sapply(ann.mod$module, function(m) assignChemTax(m)))
  ann.mod[ann.mod == "character(0)"] <- ""
  mod.assn$classes <- as.character(sapply(mod.assn$module, function(m) ann.mod[which(ann.mod$module == m),2]))
  mod.assn[mod.assn == "character(0)"] <- ""
  mod.assn$classes <- gsub(",",";",mod.assn$classes)
  mod.assn
}


