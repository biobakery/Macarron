#' Cluster metabolic features based on covarying abundances
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param dist.mat distance matrix from function MACARRoN::makeDisMat()
#' @param annot a feature annotation of choice e.g. Metabolite name/HMDB ID/METLIN ID
#' -such annotations are available for a handful of features in any MBX dataset
#' -must be a column name in feature annotation dataframe 
#' -must contain real IDs/names and not general annotations e.g. "Internal Standard", "redundant ion"
#' @param chem.info dataframe containing 3 columns: annot (V1), chemical sub-class (V2) and chemical class (V3)
#' @param min.mms min value for minimum module size (default: 10)
#' @param max.mms max value for minimum module size (default: 30)
#' @param brk.mms break size for calculating module sizes (default: 5)
#' with these module size parameters, properties of modules when minimum module size is 10,15,20,25, and 30 are shown.
#' @return mod.assn features grouped into modules based on covarying abundances
#' 
#' @examples 
#' w <- makeDisMat(mbx, ptype="diagnosis", preval=0.7)
#' mod.assn <- findModules(mbx, w, annot="HMDB.ID", chem.info=chem.info, min.mms=10, max.mms=30, brk.mms=5)
#' 
#' @export


findModules <- function(se, dist.mat, annot, chem.info, 
                         min.mms=10, max.mms=30, brk.mms=5){
  # Construct tree 
  tree <- hclust(as.dist(dist.mat), method="average")
  message("tree constructed")
  
  #Explore modules at different minimimum module sizes
  numlist <- seq(from=min.mms, to=max.mms, by=brk.mms)
  
  #measures of success
  sing <- NULL # singletons
  totc <- NULL # total modules
  pann <- NULL # % annotated modules
  hscc <- NULL # % features with the same ID that land up in the same module
  maxc <- NULL # max classes per module
  maxs <- NULL # max sub-classes per module
  perc <- NULL # 90th percentile classes per module
  pers <- NULL # 90th percentile sub-classes per module
  
  #package
  library(dynamicTreeCut)
  anno <-  as.data.frame(rowData(mbx))
  mod.assn <- as.data.frame(anno[colnames(dist.mat),annot])
  colnames(mod.assn) <- annot
  
  
  for(n in numlist){
    mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                               distM = as.matrix(dist.mat), 
                                                               deepSplit = TRUE, 
                                                               pamRespectsDendro = TRUE,
                                                               minClusterSize = n))
    
    # singletons
    sing <- rbind(sing, length(which(mod.assn$V2==0)))
    
    # total modules
    totc <- rbind(totc, max(mod.assn$V2))
    
    # % annotated modules
    m <- mod.assn[which(mod.assn$V2!="0"),]
    has.known <- NULL
    for (i in sort(unique(m$V2))){
      if(length(which(m[,2]== i & m[,1] != ""))!=0){
        has.known <- rbind(has.known,i)
      }
    }
    pann <- rbind(pann,length(has.known)*100/max(m$V2))

    # % features with the same ID that land up in the same module
    same.mod <- NULL
    diff.mod <- NULL
    times.obs <- data.frame(table(as.character(m[which(m[,1]!=""),1]))) #number of times a metabolite is observed
    times.obs <- times.obs[which(times.obs$Freq > 1),] #metabolites observed more than once
    for (t in times.obs$Var1){
      if(length(unique(m[which(m[,1]==t),2])) > 1){
        diff.mod <- rbind(diff.mod,t)}else{
          same.mod <- rbind(same.mod,t)
        }
    }
    hscc <- rbind(hscc,length(same.mod)*100/(length(same.mod)+length(diff.mod)))

    # functional heterogeneity of modules
    cls <- NULL
    scl <- NULL
    for (i in sort(unique(m$V2))){
      #ids <- as.character(m[which(m$V2==i),1])
      m2 <- m[which(m[,2]== i & m[,1] != ""),]
      info <- NULL
      for (d in sort(unique(m2[,1]))){
        c <- chem.info[which(chem.info[,1]==d),]
        info <- rbind(info,c)
      }
      scl <- rbind(scl, length(unique(info$V2)))
      cls <- rbind(cls, length(unique(info$V3)))
    }
    maxs <- rbind(maxs, max(scl))
    maxc <- rbind(maxc, max(cls))
    pers <- rbind(pers, quantile(scl,0.9))
    perc <- rbind(perc, quantile(cls,0.9))
  }
  
  #measures of success
  mos <- data.frame(cbind(totc, sing, pann, hscc, maxc, perc, maxs, pers))
  rownames(mos) <- numlist
  colnames(mos) <- c("Total modules","Singletons","% Annotated modules","% successful HMDBs","Max classes/module","90p classes/module","Max subclasses/module", "90p subclasses/module")
  write.csv(mos, file="MoS_by_MMS.csv")
  
  message("Print MOS:")
  print(mos)
  
  final.mms <- readline(prompt="Enter minimum module size of choice: ")
  mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                          distM = as.matrix(dist.mat), 
                                                          deepSplit = TRUE, 
                                                          pamRespectsDendro = TRUE,
                                                          minClusterSize = as.integer(final.mms)))
  colnames(mod.assn) <- c(annot,"module")
  rownames(mod.assn) <- colnames(dist.mat)
  mod.assn
  }
