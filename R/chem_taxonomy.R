#' Chemical taxonomy lookup table
#'
#' A dataset mapping metabolite IDs (HMDB or PubChem) to
#' chemical subclass and class, used by \code{decorateID()}.
#'
#' @format A data frame with 145837 rows and 4 variables:
#' \describe{
#'   \item{hmdb_id}{HMDB accession}
#'   \item{pubchem_compound_id}{PubChem CID (character, may be NA if only HMDB is known)}
#'   \item{sub_class}{Chemical subclass (character, may be NA)}
#'   \item{class}{Chemical class (character, may be NA)}
#' }
#' @usage data(chem_taxonomy)
"chem_taxonomy"