#' Create a SummarizedExperiment object
#'
#' @param feat_int a dataframe (features x samples) containing metabolic feature intensities (MS1 abundances)
#' i.e. metabolic feature abundance table.  
#' @param feat_anno a dataframe containing the available feature annotations (features x annotations).
#' ^^Must contain at least one annotation such as HMDB/METLIN ID/Metabolite name for 
#' the small fraction of "annotated" features. 
#' @param exp_meta a dataframe containing sample metadata (samples x metadata).
#' ^^Must contain a column "sample" which will be used for creating rownames
#' @return SummarizedExperiment object 
#' 
#' @examples 
#' mbx <- makeSumExp(feat_int,feat_anno,exp_meta)
#' 
#' @export


makeSumExp <- function(feat_int,feat_anno,exp_meta)
{
  # Create metabolic feature IDs: F[#]
  n <- as.numeric(nrow(feat_anno))
  ids <- lapply(seq_len(n), function(i) {paste0("F",i)})
  
  
  # Match sample IDs between metadata (exp_meta) and feature intensity data (feat_int).
  # Removes samples that do not have abundance values/metadata.
  rownames(exp_meta) <- exp_meta$sample
  exp_meta[exp_meta == ""] <- NA
  feat_int <- feat_int[,intersect(names(feat_int),rownames(exp_meta))]
  exp_meta <- exp_meta[intersect(names(feat_int),rownames(exp_meta)),]
  
  # Assign newly created feature IDs as rownames for intensity and annotation tables.
  rownames(feat_int) <- ids
  rownames(feat_anno) <- ids
  
  #Make SE object
  library(SummarizedExperiment)
  SummarizedExperiment::SummarizedExperiment(assays = feat_int,
                                             colData = exp_meta,
                                             rowData = feat_anno)
}

