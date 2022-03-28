#' Cluster metabolic features based on covarying abundances into modules
#' 
#' @param se SummarizedExperiment object created using Macarron::prepInput().
#' @param w distance matrix from function Macarron::makeDisMat().
#' @param input_taxonomy chemical taxonomy file with 3 columns specifying annotation, subclass and class of annotated features. 
#' Can be created using the decorateID.R utility of Macarron. 
#' Annotation specified with "standard_identifier" and annotation in the first column of the chemical taxonomy file must match.
#' @param standard_identifier name or index of column containing HMDB or PubChem IDs. Default: Column 1 in annotation dataframe.
#' @param min_module_size minimum module size to be used for module identification with dynamicTreeCut::cutreeDynamic(). 
#' Default is cube root of number of prevalent features.
#' @param evaluateMOS examine measure of success for modules identified using min_module_size, min_module_size + 5, min_module_size + 10, min_module_size - 5, min_module_size - 10
#' 
#' @return mod.assn metabolic features clustered into "modules" based on covarying abundances and measures of success.
#' 
#' @examples 
#' prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
#' abundances_df = read.csv(file = prism_abundances, row.names = 1)
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
#' metadata_df = read.csv(file = prism_metadata, row.names = 1)
#' met_taxonomy = system.file("extdata", "demo_taxonomy.csv", package="Macarron")
#' taxonomy_df = read.csv(file = met_taxonomy)
#' mbx <- Macarron::prepInput(input_abundances = abundances_df,
#'                             input_annotations = annotations_df,
#'                             input_metadata = metadata_df)
#' w <- Macarron::makeDisMat(se = mbx)
#' modules.assn <- Macarron::findMacMod(se = mbx, 
#'                                      w = w,
#'                                      input_taxonomy = taxonomy_df) 
#' 
#' 
#' 
#' @export

findMacMod <- function(se, 
                        w, 
                        input_taxonomy,
                        standard_identifier = 1,
                        min_module_size = NULL,
                        evaluateMOS = TRUE)
{
  # Construct tree 
  tree <- hclust(as.dist(w), method="average")
  message("Tree constructed")
  
  # Setting the minimum module size
  if(is.null(min_module_size)){
    min_module_size = round((nrow(w))^(1/3))
  }else{
    min_module_size = as.numeric(as.character(min_module_size))
  }
  
  # Module assignments
  anno <- as.data.frame(SummarizedExperiment::rowData(se))
  mod.assn <- as.data.frame(anno[colnames(w),standard_identifier])
  rownames(mod.assn) <- colnames(w)
  if(is.character(standard_identifier)){
    names(mod.assn) <- standard_identifier
  }else{
    names(mod.assn) <- names(anno[standard_identifier])
  }
  
  #-------------------------------------
  # Measures of success
  #-------------------------------------
  if(evaluateMOS){
    
    message("Evaluating measures of success")
    
    # Range of min_module_size for measures of success
    min_module_size.list <- c(min_module_size - 10, min_module_size - 5, min_module_size, min_module_size + 5, min_module_size + 10)
    min_module_size.list <- min_module_size.list[which(min_module_size.list > 1)]
    
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
    
    for (n in min_module_size.list){
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
          feat.tax <- input_taxonomy[which(input_taxonomy[,1] == d),]
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
                          function(x) as.character(input_taxonomy[which(input_taxonomy[,1] == x),3]))
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
    
    # Combining all measures of success
    mac.mos <- data.frame(cbind(min_module_size.list, totc, sing, pann, hscc, maxc, perc, maxs, pers, fham))
    rownames(mac.mos) <- NULL
    colnames(mac.mos) <- c("Minimum module size (min_module_size)",
                           "Total modules",
                           "Singletons",
                           "% Annotated modules",
                           "% Consistent assignments",
                           "Max classes per module",
                           "90p classes per module",
                           "Max subclasses per module",
                           "90p subclasses per module",
                           "% Features in HAM")
  }else{
    mac.mos <- "Measures of success not evaluated."
  }
  
  #-------------------------------------
  # Module assignments
  #-------------------------------------  
  mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                          distM = as.matrix(w), 
                                                          deepSplit = TRUE, 
                                                          pamRespectsDendro = TRUE,
                                                          minClusterSize = min_module_size))
  names(mod.assn)[2] <- "module"
  ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
  ann.mod <- ann.mod[which(ann.mod > 0)]
  ann.mod <- as.data.frame(ann.mod)
  colnames(ann.mod) <- "module"
  
  # Annotate modules with chemical composition information based on annotated features
  assignChemTax <- function(m){
    if(m > 0){
      annotated.features <- unique(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),1])
      classes <- toString(unique(input_taxonomy[which(input_taxonomy[,1] %in% annotated.features),3]))
    }else{
      classes <- ""
    }
  }
  ann.mod$classes <- as.character(sapply(ann.mod$module, function(m) assignChemTax(m)))
  ann.mod[ann.mod == "character(0)"] <- ""
  mod.assn$classes <- as.character(sapply(mod.assn$module, function(m) ann.mod[which(ann.mod$module == m),2]))
  mod.assn[mod.assn == "character(0)"] <- ""
  mod.assn$classes <- gsub(",",";",mod.assn$classes)
  mod.assn <- list(mod.assn,mac.mos)
  mod.assn
}


