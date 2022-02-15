#' Calculate q-value of differential abundance of metabolic features.
#' 
#' This function uses the MaAsLin2 package for estimating q-value of differential abundance.
#' Multiple fixed and random effects can be specified for fitting the multiple regression model. 
#' Default analysis method is "LM". Can be run on multiple cores.
#' metadata_variable and ref (reference group) should be the same as the one specified for effect size calculation.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp().
#' @param mod.assn the output of MACARRoN::findMacMod().
#' @param metadata_variable metadata of interest. Default: Column 2 of metadata table. 
#' Note: metadata_variable must be consistent across ava, q-value and effect-size calculations.
#' @param fixed_effects fixed effects, comma delimited e.g. c("metadata1","metadata2"). Default: all columns in metadata.
#' @param random_effects random effects, comma delimited. Default: NULL.
#' @param reference a reference level/group in each metadata, comma delimited e.g. c("metadata1,ref1","metadata2,ref2"). 
#' Default: alphabetically first phenotype/condition will be used as reference. 
#' Note: Reference must be specified for metadata with more than 2 levels.
#' @param output_folder the name of the output folder where all MaAsLin2 results will be written. Default: maaslin2_output
#' @param cores the number of R processes to run in parallel.
#' @param plot_heatmap 	Maaslin2 option-Generate a heatmap for the significant associations. Default: TRUE
#' @param plot_scatter 	Maaslin2 option-Generate scatter plots for the significant associations. Default: FALSE
#' @param heatmap_first_n Maaslin2 option-Generate heatmap for top n significant associations. Default: 50
#' @return mac.qval q-value of metabolic features in phenotypes of interest.
#' 
#' @examples
#' prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
#' abundances_df = read.csv(file = prism_abundances, row.names = 1)
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
#' metadata_df = read.csv(file = prism_metadata)
#' met_taxonomy = system.file("extdata", "demo_taxonomy.csv", package="Macarron")
#' taxonomy_df = read.csv(file = met_taxonomy)
#' mbx <- Macarron::makeSumExp(input_abundances = abundances_df,
#'                             input_annotations = annotations_df,
#'                             input_metadata = metadata_df)
#' w <- Macarron::makeDisMat(se = mbx)
#' modules.assn <- Macarron::findMacMod(se = mbx, 
#'                                      w = w,
#'                                      input_taxonomy = taxonomy_df)
#' mets.qval <- Macarron::calQval(se = mbx,
#'                                mod.assn = modules.assn)
#' 
#' @export


calQval <- function(se, 
                    mod.assn, 
                    metadata_variable = NULL,
                    fixed_effects = NULL, 
                    random_effects = NULL,
                    reference = NULL,
                    output_folder = NULL,
                    cores = 1,
                    plot_heatmap = TRUE,
                    plot_scatter = FALSE,
                    heatmap_first_n = 50
                    )
  {
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(SummarizedExperiment::assay(se))
  fint <- fint[rownames(mod.assn),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  
  # q-value estimation with Maaslin2
 
  # Default output folder
  if(is.null(output_folder)){
    output_folder = "maaslin2_output"
  }else{
    output_folder = output_folder
  }
  
  # Fitting
  fit.data <- Maaslin2::Maaslin2(input_data = fint, 
                       input_metadata = meta, 
                       output = output_folder, 
                       fixed_effects = fixed_effects, 
                       random_effects = random_effects,
                       reference = reference,
                       normalization = "NONE",
                       transform = "NONE",
                       min_prevalence = 0,
                       min_abundance = 0,
                       cores = cores,
                       plot_heatmap = plot_heatmap,
                       plot_scatter = plot_scatter,
                       heatmap_first_n = heatmap_first_n)
  
  # Adjusting p-values
  if(is.null(metadata_variable)){
    metadata_variable <- names(SummarizedExperiment::colData(se))[1]
  }else{
    metadata_variable = metadata_variable
  }
  qval <- as.data.frame(fit.data$results[which(fit.data$results$metadata == metadata_variable),
                                              c("feature","metadata","value","coef","pval")])
  mac.qval <- plyr::ddply(qval, plyr::.(metadata,value), 
                     transform, qvalue=as.numeric(stats::p.adjust(as.numeric(pval), "BH")))
  mac.qval <- mac.qval[,c("feature","metadata","value","qvalue")]
  mac.qval
}
