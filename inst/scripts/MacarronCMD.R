#!/usr/bin/env Rscript


###############################################################################

# Macarron Command Line 

###############################################################################

###############################################################################
# Set the default options
###############################################################################

execution_mode_choices <- c("serial", "multi")

args <- list()
args$input_abundances <- NULL
args$input_annotations <- NULL
args$input_metadata <- NULL
args$input_taxonomy <- NULL
args$output <- "Macarron_output"
args$metadata_variable <- 1
args$min_prevalence <- 0.7
args$execution_mode <- execution_mode_choices[1]
args$standard_identifier <- 1
args$anchor_annotation <- 2
args$min_module_size <- NULL
args$fixed_effects <- NULL
args$random_effects <- NULL
args$reference <- NULL
args$cores <- 1
args$plot_heatmap <- TRUE
args$plot_scatter <- FALSE
args$heatmap_first_n <- 50
args$show_best <- TRUE
args$priority_threshold <- 0.9
args$per_module <- 10
args$per_phenotype <- 1000
args$only_characterizable <- TRUE


###############################################################################
# Add command line arguments 
###############################################################################
options <-
  optparse::OptionParser(usage = paste(
    "%prog [Inputs]\n",
    "<feature_abundances.csv>\n",
    "<feature_annotations.csv>\n",
    "<sample_metadata.csv>\n",
    "<chemical_taxonomy.csv> "
  )
  )
options <-
  optparse::add_option(
    options,
    c("-o", "--output"),
    type = "character",
    dest = "output",
    default = args$output,
    help = paste0("Output folder name ",
                  "[Default: %default]"
    )
  )
options <- 
  optparse::add_option(
    options,
    c("-m", "--metadata_variable"),
    type = "character",
    dest = "metadata_variable",
    default = args$metadata_variable,
    help = paste("Metadata column with bioactivity-relevant epidemiological or environmental conditions (factors)",
                 "[Default: Second column of metadata CSV file]")
  )
options <-
  optparse::add_option(
    options,
    c("-p", "--min_prevalence"),
    type = "double",
    dest = "min_prevalence",
    default = args$min_prevalence,
    help = paste0("Minimum percent of samples of each condition in <metadata variable> ",
                  "in which metabolic feature is recorded ",
                  "[Default: %default]"
    )
  )
options <- 
  optparse::add_option(
    options,
    c("-e", "--execution_mode"),
    type = "character",
    dest = "execution_mode",
    default = args$execution_mode,
    help = paste0("BiocParallel Class for execution [Default: %default] [Choices:",
                  toString(execution_mode_choices),
                  "]"
    )
  )
options <- 
  optparse::add_option(
    options,
    c("-k", "--standard_identifier"),
    type = "character",
    dest = "standard_identifier",
    default = args$standard_identifier,
    help = paste("Annotation such as HMDB/METLIN/PubChem ID for an annotated metabolite feature",
                 "[Default: Second column of annotation CSV file]")
  )
options <- 
  optparse::add_option(
    options,
    c("-a", "--anchor_annotation"),
    type = "character",
    dest = "anchor_annotation",
    default = args$anchor_annotation,
    help = paste("Metabolite name for an annotated metabolite feature",
                 "[Default: Third column of annotation CSV file]")
  )
options <- 
  optparse::add_option(
    options,
    c("-c", "--min_module_size"),
    type = "double",
    dest = "min_module_size",
    default = args$min_module_size,
    help = paste("Minimum module size for module generation",
                 "[Default: cube root of number of prevalent features]")
  )
options <- 
  optparse::add_option(
    options,
    c("-f", "--fixed_effects"),
    type = "character",
    dest = "fixed_effects",
    default = args$fixed_effects,
    help = paste("Maaslin2 parameter; Fixed effects for differential abundance linear model,",
                 "comma-delimited for multiple effects",
                 "[Default: all]")
  )
options <- 
  optparse::add_option(
    options,
    c("-r", "--random_effects"),
    type = "character",
    dest = "random_effects",
    default = args$random_effects,
    help = paste("Maaslin2 parameter; Random effects for differential abundance linear model,",
                 "comma-delimited for multiple effects",
                 "[Default: none]")
  )
options <-
  optparse::add_option(
    options,
    c("-d", "--reference"),
    type = "character",
    dest = "reference",
    default = args$reference,
    help = paste("Maaslin2 parameter; The factor to use as a reference for",
                 "a variable with more than two levels",
                 "provided as a string of 'variable,reference'",
                 "comma delimited for multiple variables",
                 "[Default: Alphabetically first]"
    )
  )
options <-
  optparse::add_option(
    options,
    c("-n", "--cores"),
    type = "double",
    dest = "cores",
    default = args$cores,
    help = paste("Maaslin2 parameter; The number of R processes to",
                 "run in parallel [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-l", "--plot_heatmap"),
    type = "logical",
    dest = "plot_heatmap",
    default = args$plot_heatmap,
    help = paste("Maaslin2 parameter; Generate a heatmap for",
                 "significant associations"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-s", "--plot_scatter"),
    type = "logical",
    dest = "plot_scatter",
    default = args$plot_scatter,
    help = paste("Maaslin2 parameter; Generate scatter plots for",
                 "significant associations"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-t", "--heatmap_first_n"),
    type = "double",
    dest = "heatmap_first_n",
    default = args$heatmap_first_n,
    help = paste("Maaslin2 parameter; Generate heatmap for top n",
                 "significant associations [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-b", "--show_best"),
    type = "logical",
    dest = "show_best",
    default = args$show_best,
    help = paste("Write 1000 or fewer highly prioritized metabolic features 
                 into a separate file [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-w", "--priority_threshold"),
    type = "double",
    dest = "priority_threshold",
    default = args$priority_threshold,
    help = paste("Cut-off of priority score for showing highly
                 prioritized features [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-x", "--per_module"),
    type = "double",
    dest = "per_module",
    default = args$per_module,
    help = paste("Show first n highly prioritized features in a module [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-y", "--per_phenotype"),
    type = "double",
    dest = "per_phenotype",
    default = args$per_phenotype,
    help = paste("Show highly prioritized n features per phenotype/condition [Default: %default]"
    )
  )

options <-
  optparse::add_option(
    options,
    c("-z", "--only_characterizable"),
    type = "logical",
    dest = "only_characterizable",
    default = args$only_characterizable,
    help = paste("Show highly prioritized features in modules which contain at least
                 one annotated metabolite [Default: %default]"
    )
  )

if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  
  # check if four positional arguments are provided
  if (length(positional_args) != 4) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments"
      )
    )
  }

mac.result <- Macarron::Macarron(positional_args[1],
                                 positional_args[2],
                                 positional_args[3],
                                 positional_args[4],
                                 current_args$output,
                                 current_args$metadata_variable,
                                 current_args$min_prevalence,
                                 current_args$execution_mode,
                                 current_args$standard_identifier,
                                 current_args$anchor_annotation,
                                 current_args$min_module_size,
                                 current_args$fixed_effects,
                                 current_args$random_effects,
                                 current_args$reference,
                                 current_args$cores,
                                 current_args$plot_heatmap,
                                 current_args$plot_scatter,
                                 current_args$heatmap_first_n,
                                 current_args$show_best,
                                 current_args$priority_threshold,
                                 current_args$per_module,
                                 current_args$per_phenotype,
                                 current_args$only_characterizable)     
}
