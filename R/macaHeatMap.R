#!/usr/bin/env Rscript

###############################################################################

# MACARRoN Utility: Abundance heatmap of prioritized metabolites

###############################################################################

#load the required libraries, report an error if they are not installed

for (lib in c('grDevices', 'ComplexHeatmap', 'circlize')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

args <- list()
args$summarized_input <- NULL
args$prioritized_features <- NULL
args$top_percent <- 25
args$feature_annotation <- NULL

###############################################################################
# Add command line arguments 
###############################################################################

options <-
  optparse::OptionParser(usage = paste(
    "%prog [Inputs]",
    " <summarized_input.csv> ",
    "<prioritized_features.csv>"
    )
)
options <-
  optparse::add_option(
    options,
    c("-p", "--top_percent"),
    type = "double",
    dest = "top_percent",
    default = args$top_percent,
    help = paste0("Top-ranked p percent features to be shown ",
                  "[Default: %default]"
                  )
)
options <-
  optparse::add_option(
    options,
    c("-a", "--feature_annotation"),
    type = "character",
    dest = "feature_annotation",
    default = args$feature_annotation,
    help = paste0("Annotation to be added to features ",
                  "[Default: NULL]"
                  )
)

###############################################################################
# Main macaHeatMap function
###############################################################################

macaHeatMap <- 
  function(
    summarized_input,
    prioritized_features,
    top_percent = 25,
    feature_annotation = NULL
)
{
  # Read summarized input
  if (is.character(summarized_input)) {
      sum_input <- read.csv(summarized_input, header = TRUE, row.names = 1)
      } else {
      sum_input <- summarized_input
  }

  # Read prioritized features
  if (is.character(prioritized_features)) {
      pri_feats <- read.csv(prioritized_features, header = TRUE, row.names = 1)
      } else {
      pri_feats <- prioritized_features
  }

  # Choose top p percent features and their abundance
  #----------------------------------------------------------------------------------
  pri_feats$priority_percentage <- sapply(seq(1:nrow(pri_feats)), function(p) p*100/(nrow(pri_feats)))
  top_feats <- pri_feats[which(pri_feats$priority_percentage <= top_percent),]
  top_feats <- top_feats[order(top_feats$module),]
  only_abun <- as.data.frame(t(sum_input[,2:ncol(sum_input)]))
  colnames(only_abun) <- rownames(sum_input)
  top_abun <- as.matrix(log2(only_abun[rownames(top_feats), rownames(sum_input[order(sum_input[,1]),])] + 1))
  max_abun <- max(top_abun, na.rm = TRUE)
  min_abun <- min(top_abun, na.rm = TRUE)
  mean_abun <- mean(top_abun, na.rm = TRUE)
  color_abun <- circlize::colorRamp2(c(max_abun, mean_abun, min_abun), c("#0a0182","#d65f4f","#edca82"))

  # Annotation colors
  #----------------------------------------------------------------------------------
  # for <metadata_variable>
  metadata_colors <- c("#F16A70","#B1D877","#8CDCDA","#4D4D4D","#EDCA82","#ADD9FE")
  my_grps = unique(sum_input[,1])
  my_grps_colors = metadata_colors[sample(seq(1:6), length(unique(sum_input[,1])))]
  names(my_grps_colors) <- my_grps

  # for <prioritized_features>
  # (I) perturbation
  perturbation_colors <- c("#E36683","#94B4FF")
  names(perturbation_colors) <- sort(unique(top_feats$perturbation))
  # (II) modules
  module_colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  my_modules = unique(top_feats$module)
  my_module_colors = module_colors[sample(seq(1:433), length(unique(top_feats$module)))]
  names(my_module_colors) <- my_modules

  # Column (metadata) annotation 
  #----------------------------------------------------------------------------------
  colAnn = HeatmapAnnotation(
                           metadata = sum_input[order(sum_input[,1]),1],
                           which = 'col',
                           col = list(metadata = my_grps_colors),
                           annotation_label = names(sum_input)[1],
                           annotation_name_side = "left",
                           annotation_width = unit(c(1, 4), 'cm'),
                           gap = unit(1, 'mm')
                           )

  # Row (features) annotation 
  #----------------------------------------------------------------------------------     
  rowAnn = HeatmapAnnotation(
                           module = top_feats$module,
                           perturbation = top_feats$perturbation,
                           which = 'row',
                           col = list(module = my_module_colors, 
                                      perturbation=perturbation_colors)
                           )
  if(is.null(feature_annotation)){
    rowAnn2 = HeatmapAnnotation(
                            priority = anno_barplot(top_feats$rank, bar_width = 0.5, gp=gpar(fill = 9)),
                            annotation_label = "rank",
                            which = 'row'
                            )
  }else{                                  
    rowAnn2 = HeatmapAnnotation(
                            priority = anno_barplot(top_feats$rank, bar_width = 0.5),
                            annotation_label = "rank",
                            which = 'row',
                            std_id = anno_text(top_feats[,feature_annotation], gp = gpar(fontsize=9),)
                            )
  }
  
  # Plot Heatmap
  #---------------------------------------------------------------------------------- 
  my_heatmap = ComplexHeatmap::Heatmap(top_abun, 
        col = color_abun, 
        na_col = "white",
        name = "log2(abundance)",
        top_annotation = colAnn,
        left_annotation = rowAnn,
        right_annotation = rowAnn2,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8)
        )          

  pdf("prioritized_features_heatmap.pdf", width=8.5, height=8.5)
  ComplexHeatmap::draw(my_heatmap, annotation_legend_side = "top", heatmap_legend_side = "left")                                                   
  dev.off()
  #return(pri_heatmap)
}

if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  # check three positional arguments are provided
  if (length(positional_args) != 2) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments",
            )
    )
  }
  heatmap_object <- macaHeatMap(positional_args[1],
               positional_args[2],
               current_args$top_percent,
               current_args$feature_annotation
               )     
}