#' Rank metabolic features and prioritize based on predicted bioactivity.
#' 
#' Metabolic features are ranked based on AVA, and q-value and effect size of
#' differential abundance. The harmonic mean of these three ranks is calculated and used as 
#' the meta-rank to prioritize potentially bioactive features in a phenotype (or condition). 
#' Top-ranked features have good relative abundance, and are significantly perturbed 
#' in the specified environment/phenotype.
#' 
#' @param se SummarizedExperiment object created using MACARRoN::makeSumExp()
#' @param mod.assn the output of MACARRoN::findMacMod()
#' @param mac.ava the output of MACARRoN::calAVA()
#' @param mac.qval the output of MACARRoN::calQval()
#' @param mac.es the output of MACARRoN::calES()
#' 
#' @return mac.result - metabolic features listed according to priority 
#' 
#' @examples 
#' mac.result <- prioritize(se, mod.assn, mac.ava, mac.qval, mac.es)

prioritize <- function(se,
                       mod.assn,
                       mac.ava,
                       mac.qval,
                       mac.es){
  # packages
  suppressPackageStartupMessages(require("psych", character.only = TRUE))
  
  # Test phenotypes
  test.grps <- unique(mac.qval$value)
  
  # Prioritize for each phenotype
  prioritize.each <- function(p){
    sub.qval <- mac.qval[which(mac.qval$value == p),]
    rownames(sub.qval) <- sub.qval$feature
    all.params <- as.data.frame(cbind(rownames(mac.ava),
                                      mac.ava[,"ava"],
                                      sub.qval[rownames(mac.ava),"qvalue"],
                                      mac.es[rownames(mac.ava),p]))
    colnames(all.params) <- c("feature","ava","qval","es")
    
    
    all.params$es <- as.numeric(as.character(all.params$es))
    all.params$ava <- as.numeric(as.character(all.params$ava))
    all.params$qval <- as.numeric(as.character(all.params$qval))
    
    # Assigning direction of perturbation
    all.params$status <- ""
    all.params[which(all.params$es < 0),"status"] <- paste0("depleted in ", p)
    all.params[which(all.params$es > 0),"status"] <- paste0("enriched in ", p)
    all.params$es <- abs(all.params$es)
    
    # Ranks
    message("Assigning ranks")
    all.params$ava_rank <- rank(all.params$ava)
    all.params$qval_rank <- rank(-all.params$qval)
    all.params$es_rank <- rank(all.params$es)
    
    # Meta-rank
    message("Calculating meta-rank and prioritizing metabolic features")
    all.params$meta_rank <- harmonic.mean(t(all.params[,6:8]))
    ranked.features <- all.params[order(-all.params$meta_rank),]
    rank.perc <- ecdf(all.params$meta_rank)
    ranked.features$rank_percentile <- sapply(ranked.features$meta_rank, function(x) rank.perc(x))
    ranked.features <- as.data.frame(ranked.features)
  }
  prioritized.features <- as.data.frame(do.call(rbind, lapply(test.grps, prioritize.each)))
  prioritized.features$module <- mod.assn[prioritized.features$feature,"module"]
  prioritized.features$anchor <- mac.ava[prioritized.features$feature,"anchor"]
  prioritized.features$module_composition <- mod.assn[prioritized.features$feature,"classes"]
  prioritized.features$characterizable <- 0
  prioritized.features[which(prioritized.features$anchor != ""),"characterizable"] <- "1"
  anno <- as.data.frame(rowData(se))
  prioritized.features$annotation1 <- anno[prioritized.features$feature, 1]
  prioritized.features$annotation2 <- anno[prioritized.features$feature, 2]
  prioritized.features$ava <- round(prioritized.features$ava, 4)
  prioritized.features$es <- round(prioritized.features$es, 4)
  prioritized.features$rank_percentile <- round(prioritized.features$rank_percentile, 4)
  
  # Final table of results
  mac.result <- cbind(prioritized.features[,c("feature",
                                              "annotation1",
                                              "annotation2",
                                              "rank_percentile",
                                              "status",
                                              "module",
                                              "anchor",
                                              "module_composition",
                                              "characterizable",
                                              "ava",
                                              "qval",
                                              "es")],
                      anno[prioritized.features$feature, c(3:ncol(anno))])
  mac.result$feature <- gsub("F","",mac.result$feature)
  colnames(mac.result) <- c("Feature index",
                            names(anno)[1],
                            names(anno)[2],
                            "Rank percentile",
                            "Status",
                            "Module",
                            "Anchor",
                            "Related classes",
                            "Covaries with standard",
                            "AVA",
                            "qvalue",
                            "effect size",
                            names(anno)[3:ncol(anno)])
  write.csv(mac.result, file="prioritized_metabolites_all.csv")
  write.csv(mac.result[which(mac.result$characterizable == 1),], file="prioritized_metabolites_characterizable.csv")
  mac.result
}