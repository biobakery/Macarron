#' Calculate q-value of differential abundance of metabolic features.
#' 
#' This function uses the MaAsLin2 package for estimating q-value of differential abundance.
#' Multiple fixed and random effects can be specified for fitting the multiple regression model. 
#' Default analysis method is "LM". Can be run on multiple cores.
#' ptype and ref (reference group) should be the same as the one specified for effect size calculation.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp().
#' @param mod.assn the output of MACARRoN::findMacMod().
#' @param ptype metadata of interest. Default: Column 2 of metadata table. 
#' Note: ptype must be consistent across ava, q-value and effect-size calculations.
#' @param fix.eff fixed effects, comma delimited e.g. c("metadata1","metadata2"). Default: all columns in metadata.
#' @param ran.eff random effects, comma delimited. Default: NULL.
#' @param reference a reference level/group in each metadata, comma delimited e.g. c("metadata1,ref1","metadata2,ref2"). 
#' Default: alphabetically first phenotype/condition will be used as reference. 
#' Note: Reference must be specified for metadata with more than 2 levels.
#' @param output.folder the name of the output folder where all MaAsLin2 results will be written. Default: maaslin2_output
#' @param ncores the number of R processes to run in parallel
#' @return mac.qval
#' 
#' @examples
#' mac.qval <- calQval(se, mod.assn, ptype, fix.eff, ran.eff, reference, output.folder)
#' 
#' @import Maaslin2
#' @importFrom plyr ddply
#' @importFrom stats p.adjust
#' 
#' @export


calQval <- function(se, 
                    mod.assn, 
                    ptype = NULL,
                    fix.eff = NULL, 
                    ran.eff = NULL,
                    reference = NULL,
                    output.folder = NULL,
                    ncores = 1
                    )
  {
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(assay(se))
  fint <- fint[rownames(mod.assn),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  meta <- as.data.frame(colData(se))
  
  # q-value estimation with Maaslin2
  for (lib in c('Maaslin2','plyr','stats')) {
    requireNamespace(lib, quietly = TRUE)
  }
  
  # Default output folder
  if(is.null(output.folder)){
    output.folder = "maaslin2_output"
  }else{
    output.folder = output.folder
  }
  
  # Fitting
  fit.data <- Maaslin2::Maaslin2(input_data = fint, 
                       input_metadata = meta, 
                       output = output.folder, 
                       fixed_effects = fix.eff, 
                       random_effects = ran.eff,
                       reference = reference,
                       normalization = "NONE",
                       transform = "NONE",
                       min_prevalence = 0,
                       min_abundance = 0,
                       cores = ncores,
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
  
  # Adjusting p-values
  if(is.null(ptype)){
    ptype <- names(colData(se))[1]
  }else{
    ptype = ptype
  }
  qval <- as.data.frame(fit.data$results[which(fit.data$results$metadata == ptype),
                                              c("feature","metadata","value","coef","pval")])
  mac.qval <- plyr::ddply(qval, .(metadata,value), 
                     transform, qvalue=as.numeric(stats::p.adjust(as.numeric(pval), "BH")))
  mac.qval <- mac.qval[,c("feature","metadata","value","qvalue")]
  mac.qval
}