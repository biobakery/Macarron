#' Create a SummarizedExperiment object
#'
#' @param feat_int a dataframe (features x samples) containing metabolic feature abundances (MS1 abundances)
#' @param feat_anno a dataframe containing the available feature annotations (features x annotations).
#' ^^Column 1 must contain at least one annotation such as HMDB ID/Pubchem CID/Metabolite name for 
#' the identified/annotated features in column 1. 
#' @param exp_meta a dataframe containing sample metadata (samples x metadata).
#' ^^Column 1 will be used for creating rownames
#' ^^Column 2 should identify phenotypes or conditions (categorical metadata) associated with the samples. 
#' Must not contain NA. Rows with no specified phenotype/condition will be removed.
#' 
#' @return SummarizedExperiment object 
#' 
#' @examples 
#' se <- makeSumExp(feat_int,feat_anno,exp_meta)
#' 
#' @export


makeSumExp <- function(feat_int,feat_anno,exp_meta)
{
  # Create metabolic feature IDs: F[#]
  n <- as.numeric(nrow(feat_anno))
  ids <- lapply(seq_len(n), function(i) {paste0("F",i)})
  
  
  # Sample names
  rownames(exp_meta) <- exp_meta[,1]
  exp_meta <- exp_meta[,-1]
  exp_meta[exp_meta == ""] <- NA
  
  # Remove samples without phenotype or condition metadata
  exp_meta <- exp_meta[which(!is.na(exp_meta[1])),]
  
  # Removes samples that do not have abundance values/metadata.
  feat_int <- feat_int[,intersect(names(feat_int),rownames(exp_meta))]
  exp_meta <- exp_meta[intersect(names(feat_int),rownames(exp_meta)),]
  
  message(paste0("Samples with both abundances and metadata: ",nrow(exp_meta)))
  
  # Assign newly created feature IDs as rownames for intensity and annotation tables.
  rownames(feat_int) <- ids
  rownames(feat_anno) <- ids
  
  #Make SE object
  suppressPackageStartupMessages(require("SummarizedExperiment", character.only = TRUE))
  se <- SummarizedExperiment::SummarizedExperiment(assays = feat_int,
                                                   colData = exp_meta,
                                                   rowData = feat_anno)
  message("Summarized experiment created")
  se
}

