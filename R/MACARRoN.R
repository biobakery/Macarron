#' Macarron
#' 
#' @param input_abundances a comma-delimited file or dataframe (features x samples) containing metabolic feature intensities (abundances).
#' @param input_annotations a comma-delimited file or dataframe (features x annotations) containing available feature annotations.
#' @param input_metadata a comma-delimited file or dataframe (samples x metadata) containing sample metadata.
#' @param input_taxonomy a comma-delimited file or dataframe containing the chemical class and subclass information of annotated features.
#' @param output name of the folder where Macarron output files will be written. Default: "Macarron_output".
#' @param metadata_variable Name or index of the column that identifies the phenotypes/conditions in the study. Default: Column 1 of metadata dataframe.
#' @param min_prevalence prevalence threshold (percentage). Default = 0.7.
#' @param execution_mode BiocParallel execution mode. Options: "serial" or "multi" Default = "serial".
#' @param standard_identifier Name or index of column containing HMDB or PubChem IDs. Default: Column 1 in annotation dataframe.
#' @param anchor_annotation Name or index of column containing common names of the annotated metabolite. Default: Column 2 of annotation dataframe.
#' @param min_module_size Integer that defines the size of the smallest covariance module. Default: Cube root of number of prevalent metabolic features. 
#' @param fixed_effects Covariates for linear modeling with MaAsLin2. Default: All columns of metadata dataframe.
#' @param random_effects Random effects for linear modeling with MaAsLin2. Default: NULL.
#' @param reference Reference category (factor) in categorical metadata covariates containing three or more levels. Must be provided as a string of 'covariate,reference' semi-colon delimited for multiple covariates.
#' @param cores MaAsLin2 option-The number of R processes to be run in parallel.
#' @param plot_heatmap 	MaAslin2 option-Generate a heatmap for the significant associations. Default: TRUE
#' @param plot_scatter 	MaAslin2 option-Generate scatter plots for the significant associations. Default: FALSE
#' @param heatmap_first_n MaAslin2 option-Generate heatmap for top n significant associations. Default = 50
#' @param show_best write 1000 or fewer highly prioritized metabolic features into a separate file. Default: TRUE
#' @param priority_threshold cut-off of priority score for showing highly prioritized features. Default = 0.9
#' @param per_module show first n highly prioritized features in a module. Default = 10 
#' @param per_phenotype show highly prioritized n features per phenotype/condition. Default = 1000
#' @param only_characterizable show highly prioritized features in modules which contain at least one annotated metabolite. Default = TRUE
#' 
#' @return mac.result dataframes containing metabolic features listed according to their priority (potential bioactivity) in a phenotype of interest.
#' 
#' 
#' @examples
#' prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
#' met_taxonomy = system.file("extdata", "demo_taxonomy.csv", package="Macarron")
#' mets.prioritized <- Macarron::Macarron(input_abundances = prism_abundances,
#'                                        input_annotations = prism_annotations,
#'                                        input_metadata = prism_metadata,
#'                                        input_taxonomy = met_taxonomy)
#' 
#' @export 


Macarron <- 
  function(
    input_abundances, 
    input_annotations, 
    input_metadata, 
    input_taxonomy,
    output = "Macarron_output",
    metadata_variable = 1,
    min_prevalence = 0.7,
    execution_mode = "serial",
    standard_identifier = 1,
    anchor_annotation = 2,
    min_module_size = NULL,
    fixed_effects = NULL,
    random_effects = NULL,
    reference = NULL,
    cores = 1,
    plot_heatmap = TRUE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    show_best = TRUE,
    priority_threshold = 0.9,
    per_module = 10,
    per_phenotype = 1000,
    only_characterizable = TRUE
  )
  {
    
    # Read in the abundances, annotations, metadata and chemical classification files
    #----------------------------------------------------------------------------------
    # If a character string then this is a file name, else a data frame
    # Feature abundances
    if (is.character(input_abundances)) {
      feat_int <- data.frame(data.table::fread(input_abundances, header = TRUE, sep = ","))
      feat_int <- feat_int[,-1]
    } else {
      feat_int <- input_abundances
    }
    
    # Feature annotations
    if (is.character(input_annotations)) {
      feat_anno <- data.frame(data.table::fread(input_annotations, header = TRUE, sep = ","))
      feat_anno <- feat_anno[,-1]
    } else {
      feat_anno <- input_annotations
    }
    
    # Sample metadata
    if (is.character(input_metadata)) {
      exp_meta <- data.frame(data.table::fread(input_metadata, header = TRUE, sep = ","))
      rownames(exp_meta) <- exp_meta[,1]
      exp_meta <- exp_meta[,-1]
    } else {
      exp_meta <- input_metadata
    }
    
    # Chemical classification
    if (is.character(input_taxonomy)) {
      chem_tax <- data.frame(data.table::fread(input_taxonomy, header = TRUE, sep = ","))
    } else {
      chem_tax <- input_taxonomy
    } 
    
    
    # Create the output folder
    #----------------------------------------------------------------------------------
    if (!file.exists(output)) {
      logging::loginfo("Creating output folder.")
      dir.create(output)
    }  
    
    
    # Create log file
    #----------------------------------------------------------------------------------
    mac_log_file <- file.path(output, "Macarron.log")
    # remove log file if it exists 
    if (file.exists(mac_log_file)) {
      print(paste("Warning: Deleting existing log file: ", mac_log_file))
      unlink(mac_log_file)
    }
    logging::basicConfig(level = 'FINEST')
    logging::addHandler(logging::writeToFile, 
                        file = mac_log_file, level = "DEBUG")
    logging::setLevel(20, logging::getHandler('basic.stdout'))
    
    # Logging arguments
    #----------------------------------------------------------------------------------
    logging::loginfo("Writing function arguments to log file")
    logging::logdebug("Function arguments")
    if (is.character(input_abundances)) {
      logging::logdebug("Input abundance file: %s", input_abundances)
    }
    if (is.character(input_annotations)) {
      logging::logdebug("Input annotation file: %s", input_annotations)
    }
    if (is.character(input_metadata)) {
      logging::logdebug("Input metadata file: %s", input_metadata)
    }
    if (is.character(input_taxonomy)) {
      logging::logdebug("Input chemical classification file: %s", input_taxonomy)
    }
    if (is.character(output)) {
      logging::logdebug("Output folder: %s", output)
    }
    logging::logdebug("Phenotypic/Environmental metadata: %s", metadata_variable)
    logging::logdebug("Minimum prevalence: %f", min_prevalence)
    logging::logdebug("Execution mode: %s", execution_mode)
    logging::logdebug("Public database ID: %s", standard_identifier)
    logging::logdebug("Metabolite name: %s", anchor_annotation)
    logging::logdebug("Minimum module size: %s", min_module_size)
    logging::logdebug("Fixed effects: %s", fixed_effects)
    logging::logdebug("Random effects: %s", random_effects)
    logging::logdebug("Reference condition/level: %s", reference)
    logging::logdebug("Cores: %d", cores)
    
    
    #================================================================================== 
    # Create Summarized Experiment object (se)
    #==================================================================================
    se <- prepInput(input_abundances = feat_int,
                    input_annotations = feat_anno,
                    input_metadata = exp_meta)
    logging::loginfo("Summarized Experiment created.") 
    
    #================================================================================== 
    # Create distance matrix
    #==================================================================================
    if(is.na(as.numeric(metadata_variable))){
      metadata_variable <- metadata_variable
    }else{
      metadata_variable <- as.numeric(metadata_variable)
    }
    
    if(is.character(metadata_variable)){
      chosen_mv <- metadata_variable
    }else{
      chosen_mv <- names(SummarizedExperiment::colData(se))[metadata_variable]
    }
    logging::loginfo("Metadata chosen for prevalence filtering: %s", chosen_mv)
    
    w <- makeDisMat(se = se,
                    metadata_variable = metadata_variable,
                    min_prevalence = min_prevalence,
                    execution_mode = execution_mode)
    logging::loginfo(paste0("Distance matrix with ",nrow(w)," metabolic features created."))
    
    #================================================================================== 
    # Assign metabolic features to modules
    #================================================================================== 
    message("Initiating module detection")
    
    if(is.na(as.numeric(standard_identifier))){
      standard_identifier <- standard_identifier
    }else{
      standard_identifier <- as.numeric(standard_identifier)
    }
    
    if(is.null(min_module_size)){
      mms = round((nrow(w))^(1/3))
    }else{
      mms = as.numeric(as.character(min_module_size))
    }
    logging::loginfo("Minimum module size used for this dataset: %d", mms)
    
    mac.mod <- findMacMod(se = se,
                           w = w,
                           input_taxonomy = chem_tax,
                           standard_identifier = standard_identifier,
                           min_module_size = mms,
                           evaluateMOS = TRUE)
    modules.df <- mac.mod[[1]]
    logging::loginfo("Total number of modules detected: %d", max(modules.df[,2]))
    write.csv(mac.mod[[2]], file = file.path(output, "modules_measures_of_success.csv"))

    #================================================================================== 
    # Abundance versus anchor (AVA) calculation
    #==================================================================================   
    message("Initiating AVA calculations")
    
    if(is.na(as.numeric(anchor_annotation))){
      anchor_annotation <- anchor_annotation
    }else{
      anchor_annotation <- as.numeric(anchor_annotation)
    }
    
    mac.ava <- calAVA(se =se,
                      mod.assn = mac.mod,
                      metadata_variable = metadata_variable,
                      anchor_annotation = anchor_annotation)
    
    #================================================================================== 
    # Calculating q-value of differential abundance with MaAsLin 2
    #================================================================================== 
    message("Initiating q-value calculations")
    folder_name <- file.path(output, "maaslin2_results")
    logging::logdebug(paste0("Writing MaAsLin2 results to folder: ",folder_name))
    
    mac.qval <- calQval(se = se,
                        mod.assn = mac.mod,
                        metadata_variable = metadata_variable,
                        fixed_effects = fixed_effects, 
                        random_effects = random_effects,
                        reference = reference,
                        output_folder = folder_name,
                        cores = cores,
                        plot_heatmap = plot_heatmap,
                        plot_scatter = plot_scatter,
                        heatmap_first_n = heatmap_first_n
                        )
    
    #================================================================================== 
    # Calculating effect size of differential abundance
    #================================================================================== 
    
    message("Initiating effect size calculations")
    
    mac.es <- calES(se = se,
                    mac.qval = mac.qval)
    
    #================================================================================== 
    # Integrating ranks and prioritizing bioactives in each non-control condition
    #================================================================================== 
    message("Initiating prioritization")
    
    mac.result <- prioritize(se = se,
                             mod.assn = mac.mod,
                             mac.ava = mac.ava,
                             mac.qval = mac.qval,
                             mac.es = mac.es)
    
    # Write all prioritized metabolites
    logging::addHandler(logging::writeToFile, 
                        file = mac_log_file, level = "DEBUG")
    
    file_name <- "prioritized_metabolites_all.csv"
    file_loc = file.path(output, file_name)
    logging::loginfo(paste0("Writing all prioritized metabolites to file: ",file_loc))
    write.csv(mac.result[[1]], file=file_loc, row.names=FALSE)
    
    # Write characterizable prioritized metabolites
    file_name <- "prioritized_metabolites_characterizable.csv"
    file_loc = file.path(output, file_name)
    logging::loginfo(paste0("Writing characterizable prioritized metabolites to file: ",file_loc))
    write.csv(mac.result[[2]], file=file_loc, row.names=FALSE)
    
    #================================================================================== 
    # Highly prioritized features in each module
    #================================================================================== 
    
    if(show_best){
      best.mets <- showBest(mac.result,
                            priority_threshold = priority_threshold,
                            per_module = per_module,
                            per_phenotype = per_phenotype,
                            only_characterizable = only_characterizable)
    }
    phenotypes <- unique(best.mets$phenotype[best.mets$phenotype!=""])
    for(p in phenotypes){
      file_name <- paste0("highly_prioritized_per_module_in_",p,".csv")
      file_loc = file.path(output, file_name)
      logging::loginfo(paste0("Writing highly prioritized metabolites in ",p," to file: ",file_loc))
      write.csv(best.mets[which(best.mets$phenotype == p),], file=file_loc, row.names=FALSE)
    }
    mac.result                                                                                                                                                                                                
  }
