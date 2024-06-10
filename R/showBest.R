#' View highly prioritized bioactives grouped by modules. 
#'
#' Modules are listed in the order of priority. Only the top-ranked n features in each module 
#' are shown. The priority of a module is the ratio of number of features in it that are ranked higher than
#' the cut-off and the size of the module. This utility function makes it easier to understand default
#' prioritization results of large datasets where a few hundred metabolic features are highly-prioritized.
#' 
#' @param mac.result the output of Macarron::Macarron() or Macarron::prioritize().
#' @param priority_threshold cut-off of priority score. Default = 0.9.
#' @param per_module show first n highly prioritized features in a module. Default = 10  
#' @param per_phenotype show highly prioritized n features per phenotype/condition. Default = 1000
#' @param only_characterizable show highly prioritized features in modules which contain at least one annotated metabolite. Default = TRUE
#' 
#' @return best.mets -highly-prioritized bioactives in each module in each phenotype
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
#' mets.qval <- Macarron::calQval(se = mbx,
#'                                mod.assn = modules.assn)
#' mets.es <- Macarron::calES(se = mbx,
#'                            mac.qval = mets.qval)
#' mets.prioritized <- Macarron::prioritize(se = mbx,
#'                                          mod.assn = modules.assn,
#'                                          mac.ava = mets.ava,
#'                                          mac.qval = mets.qval,
#'                                          mac.es = mets.es)  
#' best.mets <- Macarron::showBest(mac.result = mets.prioritized)
#'                                                                 
#' 
#' @export

showBest <- function(mac.result,
                   priority_threshold = 0.9,
                   per_module = 10,
                   per_phenotype = 1000,
                   only_characterizable = TRUE){
  
  # Only modules
  mac.result <- mac.result[[1]]
  mac.result <- mac.result[which(mac.result$Module != 0),]
  mac.result$phenotype <- gsub(".*in ","",mac.result$Status)
  module_size <- as.data.frame(table(mac.result$Module))
  module_size[,2] <- module_size[,2]/2
  names(module_size) <- c("module","size")
  # Calculating delta property
  if(is.numeric(mac.result[,4])){
    anchor_values <- unique(mac.result[which(mac.result$Metabolite == mac.result$Anchor & mac.result$AVA == 1),c("Module", names(mac.result)[4])])
    rownames(anchor_values) <- anchor_values$Module
    anchor_values$Module <- as.character(anchor_values$Module)
    mac.result$delta_value <- anchor_values[as.character(mac.result$Module), 2]
    mac.result$delta_value <- mac.result[,4] - mac.result$delta_value
    mac.result$delta_value <- round(mac.result$delta_value, 3)
    delta_name <- paste0(names(mac.result)[4],"_vs_Anchor")
    names(mac.result)[ncol(mac.result)] <- delta_name
  }
  # Test phenotypes
  test.phenotypes <- unique(mac.result$phenotype)
  # Digest each phenotype
  digest.each <- function(p){
    # Subset characterizable metabolites by phenotype and priority score threshold 
    phenotype_df <- mac.result[which(mac.result$phenotype == p & 
                                 mac.result$Priority_score >= priority_threshold),]
    if(only_characterizable){
      phenotype_df <-  phenotype_df[phenotype_df$Covaries_with_standard == 1,]
    }
    module_rep <- as.data.frame(table(phenotype_df$Module))
    names(module_rep) <- c("module","primets")
    module_rep <- merge(module_rep, module_size, by="module")
    module_rep$score <- module_rep$primets/module_rep$size
    module_rep <- module_rep[order(-module_rep$score),]
    # Function to find top-ranked features and standards in each module
    analyze.module <- function(m){
        df_sub <- phenotype_df[which(phenotype_df$Module == m),]
      if(nrow(df_sub) >= per_module){
        df_pri <- head(df_sub, per_module)
      }else{
        df_pri <- df_sub
      }
      df_all <- rbind(mac.result[which(mac.result[,2] != "" & mac.result$Module == m & mac.result$phenotype == p),],df_pri, NA)
      df_all[!duplicated(df_all),]
    }
    # Top-ranked features in each phenotype
    df_cat <- as.data.frame(do.call(rbind, lapply(module_rep$module, analyze.module)))
    rownames(df_cat) <- seq_len(nrow(df_cat))
    df_cat <- head(df_cat, as.numeric(tail(rownames(head(df_cat[complete.cases(df_cat),],per_phenotype)),1)))
    df_cat <- df_cat[,c(1:4, ncol(df_cat), 5:(ncol(df_cat)-1))]
    df_cat <- cbind(seq_len(nrow(df_cat)), df_cat)
    names(df_cat)[1] <- "Row"
    df_cat[is.na(df_cat)] <- ""
    df_cat$phenotype <- p
    df_cat
  }
  best.mets <- as.data.frame(do.call(rbind, lapply(test.phenotypes, digest.each))) 
  best.mets
}