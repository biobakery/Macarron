#' Create a SummarizedExperiment object
#'
#' @param feat_int a dataframe (features x samples) containing metabolic feature intensities (abundances).
#' @param feat_anno a dataframe (features x annotations) containing the available feature annotations.
#' ^^Column 1 must contain standard annotations such as HMDB ID or Pubchem CID for 
#' the subset of identified/annotated features. 
#' @param exp_meta a dataframe (samples x metadata) containing sample metadata.
#' ^^Column 1 must identify samples. It will be used for creating rownames.
#' ^^Column 2 must identify phenotypes or conditions (categorical metadata) associated with the samples. 
#' Must not contain NA. Rows with no specified phenotype/condition will be removed.
#' 
#' @return SummarizedExperiment object 
#' 
#' @examples 
#' se <- makeSumExp(feat_int,feat_anno,exp_meta)
#' 
#' @import SummarizedExperiment
#' 
#' @export


makeSumExp <- function(feat_int,feat_anno,exp_meta)
{
  # Create metabolic feature IDs: F[#]
  n <- as.numeric(nrow(feat_anno))
  ids <- lapply(seq_len(n), function(i) {paste0("F",i)})
  
  # Assign newly created feature IDs as rownames for intensity and annotation tables.
  rownames(feat_int) <- ids
  rownames(feat_anno) <- ids
  
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
  
  #Make SE object
  se <- SummarizedExperiment::SummarizedExperiment(assays = feat_int,
                                                   colData = exp_meta,
                                                   rowData = feat_anno)
  message("Summarized experiment created")
  se
}

