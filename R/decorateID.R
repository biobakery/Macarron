#' Create a chemical taxonomy table for annotated metabolic features.
#' Uses a pregenerated chemical taxonomy table shipped with the package
#' to assign chemical subclass and class to annotated features.
#' 
#' @param input_annotations a dataframe (features x annotations) containing the available feature annotations.
#' ^^Column 1 must contain standard annotations such as HMDB ID or PubChem CID for 
#' the subset of identified/annotated metabolic features. 
#' 
#' @param chemical_taxonomy Optional data frame with columns
#'  \code{hmdb_id}, \code{pubchem_compound_id}, \code{sub_class},
#'  and \code{class}. If \code{NULL} (default), the built-in
#'   \code{chem_taxonomy} dataset is loaded via \code{data()}.
#'   
#' @return A data frame with three columns: ID (same type as the first
#'  column of \code{input_annotations}), \code{Sub_Class} and \code{Class},
#'  for the subset of features whose IDs are present in
#'  \code{chemical_taxonomy}.
#' 
#' @examples 
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' input_taxonomy <- decorateID(annotations_df)
#' 
#' @export


decorateID <- function(input_annotations,
                       chemical_taxonomy = NULL)
{
  # load internal taxonomy if user didn't supply one
  if (is.null(chemical_taxonomy)) {
    data("chem_taxonomy", package = "Macarron", envir = environment())
    chemical_taxonomy <- chem_taxonomy
  }
  
  # unique, non-empty IDs from first annotation column
  ID_list <- unique(input_annotations[, 1])
  ID_list <- ID_list[ID_list != "" & !is.na(ID_list)]
  
  # decide ID type
  is_hmdb <- all(grepl("^HMDB", ID_list))
  
  if (is_hmdb) {
    tax_df <- chemical_taxonomy[
      chemical_taxonomy$hmdb_id %in% ID_list,
      c("hmdb_id", "sub_class", "class"),
      drop = FALSE
    ]
  } else {
    tax_df <- chemical_taxonomy[
      chemical_taxonomy$pubchem_compound_id %in% ID_list,
      c("pubchem_compound_id", "sub_class", "class"),
      drop = FALSE
    ]
  }
  
  names(tax_df) <- c(names(input_annotations)[1], "Sub_Class", "Class")
  rownames(tax_df) <- NULL
  as.data.frame(tax_df)
}
