#' Calculate q-value of differential abundance of metabolic features.
#' 
#' This function uses the MaAsLin2 package for estimating q-value of differential abundance.
#' Multiple fixed and random effects can be specified for fitting the multiple regression model. 
#' Default analysis method is "LM". Can be run on multiple cores.
#' ptype and ref (reference group) should be the same as the one specified for effect size calculation.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param mod.assn the output of MACARRoN::findModules(); a dataframe containing feature names, 
#' associate annotations and module assignments
#' @param ptype a reference meta.data of choice. 
#' @param ref a reference level/group in ptype
#' @param fix.eff fixed effects, comma delimited i.e. c("a", "b")
#' @param ran.eff random effects, comma delimited
#' @param folder.name the name of the output folder where all MaAsLin2 results will be written
#' @param ncores the number of R processes to run in parallel
#' @return my.qvalue
#' 
#' @examples
#' my.effsiz <- calEffSize(mbx, mod.assn, ptype, ref)
#' my.qvalue <- calQvalue(mbx, mod.assn, ptype, ref, fix.eff, ran.eff, folder.name, ncores)
#' 
#' @export


calQvalue <- function(se, mod.assn, ptype, ref, fix.eff, ran.eff, folder.name, ncores){
  mod.assn <- as.data.frame(mod.assn)
  fint <- as.data.frame(assay(se))
  fint <- fint[rownames(mod.assn),]
  fint[is.na(fint)] <- 0
  fint <- log2(fint + 1)
  meta <- as.data.frame(colData(se))
  meta[,ptype] <- relevel(as.factor(meta[,ptype]), ref=ref)
  
  # Mpdel fitting with MaAsLin2
  library(Maaslin2)
  fit.data <- Maaslin2(input_data = fint, 
                       input_metadata = meta, 
                       output = folder.name, 
                       fixed_effects = fix.eff, 
                       random_effects = ran.eff,
                       normalization = "NONE",
                       transform = "NONE",
                       min_prevalence = 0,
                       min_abundance = 0,
                       cores = ncores,
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
  
  library(plyr)
  my.qvalue <- as.data.frame(fit.data$results[which(fit.data$results$metadata == ptype),
                                              c("feature","metadata","value","coef","pval")])
  my.qvalue <- ddply(my.qvalue, .(metadata,value), 
                     transform, qvalue=as.numeric(p.adjust(as.numeric(pval), "BH")))
  my.qvalue
}


