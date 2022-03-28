#' Calculate effect size of differential abundance of metabolic features.
#' 
#' Effect size of a metabolic feature is the difference in mean log2 transformed abundances in 
#' test and control (reference) samples. For the specified metadata variable, effect size 
#' is calculated for all test categories against the reference category.
#' 
#' @param se SummarizedExperiment object created using Macarron::prepInput().
#' @param mac.qval the output of Macarron::calQval().
#' @return mac.es effect sizes of metabolic features in phenotypes of interest.
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
#' mets.qval <- Macarron::calQval(se = mbx,
#'                                mod.assn = modules.assn)
#' mets.es <- Macarron::calES(se = mbx,
#'                            mac.qval = mets.qval)                               
#' 
#' @export

calES <- function(se,
                  mac.qval)
{
  # Abundance matrix
  fint <- as.data.frame(SummarizedExperiment::assay(se))
  fint <- fint[unique(mac.qval$feature),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  fint <- t(fint)
  
  # Get phenotype from mac.qval
  metadata_variable <- unique(mac.qval$metadata)
  phenotypes <- unique(se[[metadata_variable]])
  
  # Mean abundance of feature in each group of metadata_variable
  get.mean <- function(g){
    message(paste0("Calculating mean abundance in: ",g))
    ind <- se[[metadata_variable]] == g
    m <- sapply(colnames(fint), function(f) mean(fint[ind,f]))
  }
  all.means <- as.data.frame(do.call(cbind, lapply(phenotypes, get.mean)))
  colnames(all.means) <- phenotypes
  
  # Effect size
  test.phenotypes <- unique(mac.qval$value)
  ref.phenotype <- setdiff(phenotypes, test.phenotypes)
  get.es <- function(test.phenotype){
    message(paste0("Calculating effect size in: ",test.phenotype))
    es <- all.means[,test.phenotype]- all.means[,ref.phenotype]
  }
  mac.es <- as.data.frame(do.call(cbind, lapply(test.phenotypes, get.es)))
  colnames(mac.es) <- test.phenotypes
  rownames(mac.es) <- colnames(fint)
  mac.es
}
