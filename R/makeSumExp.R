#' Create a SummarizedExperiment object
#'
#' @param input_abundances a dataframe (features x samples) containing metabolic feature intensities (abundances).
#' @param input_annotations a dataframe (features x annotations) containing the available feature annotations.
#' ^^Column 1 must contain standard annotations such as HMDB ID or Pubchem CID for 
#' the subset of identified/annotated features. 
#' ^^Column 2 must contain metabolite name.
#' @param input_metadata a dataframe (samples x metadata) containing sample metadata.
#' ^^Column 1 must identify samples. It will be used for creating rownames.
#' ^^Column 2 must identify phenotypes or conditions (categorical metadata) associated with the samples. 
#' Must not contain NA. Rows with no specified phenotype/condition will be removed.
#' 
#' @return SummarizedExperiment object 
#' 
#' @examples
#' prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
#' abundances_df = read.csv(file = prism_abundances, row.names = 1)
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
#' metadata_df = read.csv(file = prism_metadata)
#' mbx <- Macarron::makeSumExp(input_abundances = abundances_df,
#'                             input_annotations = annotations_df,
#'                             input_metadata = metadata_df)
#' 
#' @import SummarizedExperiment
#' 
#' @export


makeSumExp <- function(input_abundances,input_annotations,input_metadata)
{
  # Create metabolic feature IDs: F[#]
  n <- as.numeric(nrow(input_annotations))
  ids <- lapply(seq_len(n), function(i) {paste0("F",i)})
  
  # Assign newly created feature IDs as rownames for intensity and annotation tables.
  rownames(input_abundances) <- ids
  rownames(input_annotations) <- ids
  
  # Sample names
  rownames(input_metadata) <- input_metadata[,1]
  input_metadata <- input_metadata[,-1]
  input_metadata[input_metadata == ""] <- NA
  
  # Remove samples without phenotype or condition metadata
  input_metadata <- input_metadata[which(!is.na(input_metadata[1])),]
  
  # Removes samples that do not have abundance values/metadata.
  input_abundances <- input_abundances[,intersect(names(input_abundances),rownames(input_metadata))]
  input_metadata <- input_metadata[intersect(names(input_abundances),rownames(input_metadata)),]
  
  message(paste0("Samples with both abundances and metadata: ",nrow(input_metadata)))
  
  #Make SE object
  se <- SummarizedExperiment::SummarizedExperiment(assays = input_abundances,
                                                   colData = input_metadata,
                                                   rowData = input_annotations)
  message("Summarized experiment created")
  se
}

