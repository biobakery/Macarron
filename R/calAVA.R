#' Calculate abundance versus anchor (AVA) of metabolic features.
#' 
#' AVA of a feature is the ratio of its abundance and the most abundant metabolite in the same module 
#' i.e. the "anchor". Anchor is an annotated/known feature if available or just the most abundant metabolic
#' feature. For every feature, mean abundance in each phenotype or condition is calculated and the maximum is 
#' considered for AVA calculation. Singletons are assigned an AVA of 1.
#' 
#' @param se SummarizedExperiment object created using Macarron::prepInput().
#' @param mod.assn the output of Macarron::findMacMod().
#' @param metadata_variable name or index of metadata column identifying phenotypes/conditions to be used for evaluating AVA. Default: Column 1 of metadata dataframe. 
#' Note: metadata_variable must be consistent across distance matrix, ava, q-value and effect-size calculations.
#' @param anchor_annotation name or index of column containing common names of the annotated metabolite. Default: Column 2 of annotation dataframe.
#' @return mac.ava abundance versus anchor values of metabolic features
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
#' mets.ava <- Macarron::calAVA(se = mbx,
#'                              mod.assn = modules.assn)
#' 
#' @export

calAVA <- function(se, 
                   mod.assn,
                   metadata_variable = 1,
                   anchor_annotation = 2)
  {
  mod.assn <- as.data.frame(mod.assn[[1]])
  fint <- as.data.frame(SummarizedExperiment::assay(se))
  fint <- fint[rownames(mod.assn),]
  fint <- t(fint)
  
  # Setting the metadata
  if(is.character(metadata_variable)){
    metadata_variable <- metadata_variable
  }else{
    metadata_variable <- names(SummarizedExperiment::colData(se))[metadata_variable]
  }
  phenotypes <- unique(se[[metadata_variable]])
  
  # Identification of the anchor in each module
  anchors <- NULL
  modules <- sort(unique(mod.assn$module))
  
  # Function for anchor and the reference phenotype
  message("Finding anchors")
  find.anchor <- function(m){
    if(dim(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])[1] > 0){
      # known features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])
    }else{
      # all features in the module
      r <- rownames(mod.assn[which(mod.assn$module == m),])
    }
    # Mean abundances in phenotypes
    means <- NULL
    for (g in phenotypes){
      ind <- se[[metadata_variable]] == g
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
    tail(means[order(means[,"mma"]),],1)
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
  message("Calculating AVA")
  cal.ava <- function(i){
    getMean <- function(g){
      ind <- se[[metadata_variable]] == g
      gmean <- mean(fint[ind,i], na.rm = TRUE)
      gmean
    }
    all.means <- vapply(phenotypes, getMean, numeric(1))
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
  anno <- as.data.frame(SummarizedExperiment::rowData(se))
  if(is.character(anchor_annotation)){
    anchor_annotation <- anchor_annotation
  }else{
    anchor_annotation <- names(anno[anchor_annotation])
  }
  
  ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
  ann.mod <- ann.mod[which(ann.mod > 0)]
  ann.mod <- as.data.frame(ann.mod)
  colnames(ann.mod) <- "module"
  
  assignAnchor <- function(m){
      anchor.feature <- rownames(ava[which(ava$module == m & ava$ava == 1),])
      anchor.name <- as.character(anno[anchor.feature, anchor_annotation])
    }
  ann.mod$anchor <- as.character(sapply(ann.mod$module, function(m) assignAnchor(m)))
  ann.mod[ann.mod == "character(0)"] <- ""
  ava$anchor <- as.character(sapply(ava$module, function(m) ann.mod[which(ann.mod$module == m),2]))
  ava[ava == "character(0)"] <- ""
  ava
}
