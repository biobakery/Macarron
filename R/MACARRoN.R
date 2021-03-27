#!/usr/bin/env Rscript --vanilla

###############################################################################

# MACARRoN

###############################################################################


#load the required libraries, report an error if they are not installed

for (lib in c('optparse', 'logging', 'data.table', 'SummarizedExperiment')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

###############################################################################
# Set the default options
###############################################################################

execution_mode_choices <- c("serial", "multi")

args <- list()
args$input_abundances <- NULL
args$input_annotations <- NULL
args$input_metadata <- NULL
args$input_classification <- NULL
args$output <- "MACARRoN_output"
args$metadata_variable <- NULL
args$control_condition <- NULL
args$min_prevalence <- 0.7
args$execution_mode <- execution_mode_choices[1]
args$standard_identifier <- NULL
args$min_module_size <- NULL
args$fixed_effects <- NULL
args$random_effects <- NULL
args$reference <- NULL
args$cores <- 1


###############################################################################
# Add command line arguments 
###############################################################################

options <-
  optparse::OptionParser(usage = paste(
    "%prog [Inputs]\n",
    "<feature_abundances.csv>\n",
    "<feature_annotations.csv>\n",
    "<sample_metadata.csv>\n",
    "<chemical_classification.csv> "
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
    c("-c", "--metadata_variable"),
    type = "character",
    dest = "metadata_variable",
    default = args$metadata_variable,
    help = paste("Metadata column with bioactivity-relevant epidemiological or environmental conditions (factors)",
            "[Default: First metadata column]")
)
options <- 
  optparse::add_option(
    options,
    c("-a", "--control_condition"),
    type = "character",
    dest = "control_condition",
    default = args$control_condition,
    help = paste("Condition within <metadata variable> to be used as the reference (control)",
                 "for differential abundance analysis and prioritization",
                 "[Default: Alphabetically first category in <metadata variable>]")
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
    help = paste("Annotation such as HMDB/METLIN/MoNA accession or metabolite name for an annotated metabolite feature",
            "[Default: First column in annotation table]")
)
options <- 
  optparse::add_option(
    options,
    c("-m", "--min_module_size"),
    type = "character",
    dest = "min_module_size",
    default = args$min_module_size,
    help = paste("Minimum module size for module generation",
            "[Default: %default (will be calculated by analyzing measures of success)]")
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
            "semi-colon delimited for multiple variables [Default: NA]"
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




###############################################################################
# Main MACARRoN function (defaults same command line) 
###############################################################################

MACARRoN <- 
  function(
    input_abundances, 
    input_annotations, 
    input_metadata, 
    input_classification,
    output = "MACARRoN_output",
    metadata_variable = NULL,
    control_condition = NULL,
    min_prevalence = 0.7,
    execution_mode = "serial",
    standard_identifier = NULL,
    min_module_size = NULL,
    fixed_effects = NULL,
    random_effects = NULL,
    reference = NULL,
    cores = 1
    )
  {

    # Read in the abundances, annotations, metadata and chemical classification files
    #----------------------------------------------------------------------------------
    # If a character string then this is a file name, else a data frame
    # Feature abundances
    if (is.character(input_abundances)) {
      feat_int <- data.frame(data.table::fread(input_abundances, header = TRUE, sep = ","))
      } else {
      feat_int <- input_abundances
    }

    # Feature annotations
    if (is.character(input_annotations)) {
      feat_anno <- data.frame(data.table::fread(input_annotations, header = TRUE, sep = ","))
      } else {
      feat_anno <- input_annotations
    }

    # Sample metadata
    if (is.character(input_metadata)) {
      exp_meta <- data.frame(data.table::fread(input_metadata, header = TRUE, sep = ","))
      } else {
      exp_meta <- input_metadata
    }

    # Chemical classification
    if (is.character(input_classification)) {
      chem_info <- data.frame(data.table::fread(input_classification, header = TRUE, sep = ","))
      } else {
      chem_info <- input_classification
    } 
    
  
    # Create the output folder
    #----------------------------------------------------------------------------------
    if (!file.exists(output)) {
      logging::loginfo("Creating output folder.")
      dir.create(output)
    }  

   
    # Create log file
    #----------------------------------------------------------------------------------
    log_file <- file.path(output, "MACARRoN.log")
    # remove log file if it exists 
    if (file.exists(log_file)) {
      print(paste("Warning: Deleting existing log file: ", log_file))
      unlink(log_file)
    }
    logging::basicConfig(level = 'FINEST')
    logging::addHandler(logging::writeToFile, 
      file = log_file, level = "DEBUG")
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
    if (is.character(input_classification)) {
      logging::logdebug("Input chemical classification file: %s", input_classification)
    }
    if (is.character(output)) {
      logging::logdebug("Output folder: %s", output)
    }
    logging::logdebug("Epidemiological/Environmental metadata: %s", metadata_variable)
    logging::logdebug("Control category in chosen metadata: %s", control_condition)
    logging::logdebug("Minimum prevalence: %f", min_prevalence)
    logging::logdebug("Execution mode: %s", execution_mode)
    logging::logdebug("Public database ID or metabolite name: %s", standard_identifier)
    logging::logdebug("Minimum module size: %s", min_module_size)
    logging::logdebug("Fixed effects: %s", fixed_effects)
    logging::logdebug("Random effects: %s", random_effects)
    logging::logdebug("Reference condition/level: %s", reference)
    logging::logdebug("Cores: %d", cores)



    # Create Summarized Experiment object
    #----------------------------------------------------------------------------------
    # Create metabolic feature IDs: F[#]
    ids <- lapply(seq_len(as.numeric(nrow(feat_anno))), function(i) {paste0("F",i)})
    # Assign newly created feature IDs as rownames for intensity and annotation tables.
    rownames(feat_int) <- ids
    rownames(feat_anno) <- ids
  
    # Match sample IDs between metadata (exp_meta) and abundance data (feat_int).
    # Removes samples that do not have features/metadata.
    rownames(exp_meta) <- exp_meta$sample
    exp_meta[exp_meta == ""] <- NA
    feat_int <- feat_int[,intersect(names(feat_int),rownames(exp_meta))]
    exp_meta <- exp_meta[intersect(names(feat_int),rownames(exp_meta)),]

    logging::loginfo("Samples with both metabolic features and metadata: %d", nrow(exp_meta))
    logging::loginfo("Total number of metabolic features: %d", nrow(feat_int))

    #Summarized Experiment object
    se <- SummarizedExperiment::SummarizedExperiment(assays = feat_int,
                                                     colData = exp_meta,
                                                     rowData = feat_anno)
    logging::loginfo("Summarized Experiment created.") 


    # Filtering features based on prevalence (default = 0.7)
    #----------------------------------------------------------------------------------
    for (lib in c('WGCNA', 'DelayedArray', 'BiocParallel', 'ff', 'ffbase')) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }

    # Abundance matrix
    mat <- DelayedArray::DelayedArray(SummarizedExperiment::assay(se))

    # Phenotype i.e. groups/conditions
    if(is.null(metadata_variable)){
      ptype <- head(colnames(as.data.frame(colData(se))),1)
    }else{
      ptype <- metadata_variable
    }
    grps <- unique(se[[ptype]])

    # Function: Get features that satisfy prevalence threshold in each condition
    .getIds <- function(g){
      ind <- se[[ptype]] == g  
      smat <- mat[,ind]
      ind <- rowMeans(is.na(smat)) <= 1 - min_prevalence
      ind
    }


    # Apply on all groups/conditions
    ind <- vapply(grps, .getIds, logical(nrow(mat)))

    # Union of features
    ind <- apply(ind,1,any)
    mat <- mat[ind,]

    logging::loginfo("Metabolic features that passed prevalence threshold: %d", nrow(mat))


    #================================================================================== 
    # Assigning filtered features to modules
    #==================================================================================
    print("Initiating module assignments")

    # Create bicor-based distance matrix
    #----------------------------------------------------------------------------------
    # Compute correlation and distance matrices
    #optimize.for = c("runtime", "memory")
    #optimize.for <- match.arg(optimize.for)
    #opt.mem <- optimize.for == "memory"
    
    
    if(nrow(mat) <= 30000){
    #Function: Correlation matrix
    .getCorMat <- function(g)
    {
      inx <- se[[ptype]] == g
      gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=min_prevalence),inx]
      tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
      rownames(tmp) <- rownames(mat)
      colnames(tmp) <- colnames(gmat)
      tmp <- DelayedArray::DelayedArray(tmp)
      tmp[rownames(gmat), colnames(gmat)] <- gmat
      tmp[is.na(tmp)] <- 0
      tmp <- log2(tmp + 1)
      options(warn= -1) 
      cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs")
      options(warn = 0)
      cmat[is.na(cmat)] <- 0
      cmat[cmat < 0] <- 0
      cmat <- as(cmat, "dsyMatrix")
      cmat
    }
    if(execution_mode == "serial"){
      exe.choice <- BiocParallel::SerialParam()
    }else if(execution_mode == "multi"){
      exe.choice <- BiocParallel::MultiParam()
    }

    # Apply on all groups/conditions
    cmats <- BiocParallel::bplapply(grps, .getCorMat, BPPARAM = exe.choice)
  
    # Keep the best observed positive correlation for each pair of features
    mmat <- do.call(pmax, c(cmats, na.rm=TRUE))
    }
    else
    # Large datasets
    {for (g in grps)
    {
      message(g)
      inx <- se[[ptype]] == g
      gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=min_prevalence),inx]
      tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
      rownames(tmp) <- rownames(mat)
      colnames(tmp) <- colnames(gmat)
      tmp <- DelayedArray::DelayedArray(tmp)
      tmp[rownames(gmat), colnames(gmat)] <- gmat
      tmp[is.na(tmp)] <- 0
      tmp <- log2(tmp + 1)
      options(warn = -1)
      cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs")
      options(warn = 0)
      assign(paste0(g,"_ff"),as.ff(cmat))
      rm(cmat)
    }
    # Keep the best observed positive correlation for each pair of features
    mmat <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
    for(j in ls(pattern="_ff")){
      for(k in 1:nrow(mmat)){
        mmat[,k] <- pmax(mmat[,k],get(j)[,k],na.rm=TRUE)
      }
    }
    rownames(mmat) <- rownames(mat)
    colnames(mmat) <- rownames(mat)
    } 
    # Beta scaling to ensure power law distribution
    mmat = mmat^3
    # Distance matrix
    w = 1 - mmat
    logging::loginfo(paste0("Distance matrix with ",nrow(w)," metabolic features created."))


    # Tree construction
    #----------------------------------------------------------------------------------
    for (lib in c('dynamicTreeCut')) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }

    # Construct tree 
    tree <- hclust(as.dist(w), method="average")
    logging::loginfo("Tree constructed.")


    # Calculating measures of success
    #----------------------------------------------------------------------------------
    # Candidates for minimum module sizes
    numlist <- seq(from=5, to=30, by=5)
    # Phenotype i.e. groups/conditions
    if(is.null(standard_identifier)){
      standard_identifier <- head(colnames(as.data.frame(rowData(se))),1)
    }else{
      standard_identifier <- standard_identifier
    }

    # measures of success
    sing <- NULL # singletons
    totc <- NULL # total modules
    pann <- NULL # % annotated modules
    fann <- NULL # % features associated with a standard 
    hscc <- NULL # % features with the same ID that land up in the same module
    maxc <- NULL # max classes per module
    maxs <- NULL # max sub-classes per module
    perc <- NULL # 90th percentile classes per module
    pers <- NULL # 90th percentile sub-classes per module

    anno <-  as.data.frame(rowData(se)) # feature annotations
    mod.assn <- as.data.frame(anno[colnames(w),standard_identifier])
    colnames(mod.assn) <- standard_identifier
    rownames(mod.assn) <- colnames(w)

    mod.assn$module <- 0         

    for(n in numlist){
      logging::loginfo(paste0("Calculating measures of success for Minimum Module Size (MMS): ",n))
      mod.assn$module <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                               distM = as.matrix(w), 
                                                               deepSplit = TRUE, 
                                                               pamRespectsDendro = TRUE,
                                                               minClusterSize = n,
                                                               verbose = 0)) 

      # singletons
      sing <- rbind(sing, length(which(mod.assn$module == 0))) 

      #total modules
      totc <- rbind(totc, max(mod.assn$module))

      # % annotated modules
      has.std <-  sort(unique(mod.assn[which(mod.assn[,1] != "" & mod.assn$module !=0 ),2]))
      pann <- rbind(pann, round(length(has.std)*100/max(mod.assn$module),2))

      # % features associated with a standard (characterizable features)
      fann <- rbind(fann, round((nrow(mod.assn[which(mod.assn[,1] == "" & mod.assn$module %in% has.std),]))*100/nrow(mod.assn),2)) 

      # %  features with the same ID that land up in the same module
      same.mod <- NULL
      diff.mod <- NULL
      std.occ <- data.frame(table(mod.assn[which(mod.assn[,1] != ""),1]))
      std.multi <- std.occ[which(std.occ$Freq > 1),1]
      for (m in std.multi){
        if(length(unique(mod.assn[which(mod.assn[,1]==m),2])) > 1){
          diff.mod <- rbind(diff.mod,m)
        }else{
          same.mod <- rbind(same.mod,m)
        }
      }
      hscc <- rbind(hscc, round(length(same.mod)*100/(length(same.mod)+length(diff.mod)),2))

      #functional heterogeneity of modules
      chem.cla <- NULL
      chem.scl <- NULL
      for (i in seq(1:max(mod.assn$module))){
        ids <- mod.assn[which(mod.assn$module == i & mod.assn[,1] != ""),1]
        ids_info <- NULL
        for (j in ids){
          c <- chem_info[which(chem_info[,1]==j),]
          ids_info <- rbind(ids_info, c)
        }
        chem.cla <- rbind(chem.cla, length(unique(ids_info$Class)))
        chem.scl <- rbind(chem.scl, length(unique(ids_info$Sub_Class)))
      }
       maxs <- rbind(maxs, max(chem.scl))
       maxc <- rbind(maxc, max(chem.cla))
       pers <- rbind(pers, quantile(chem.scl,0.9))
       perc <- rbind(perc, quantile(chem.cla,0.9))
    }

    # Compile measures of success
    mos <- data.frame(cbind(totc, sing, pann, fann, hscc, maxc, perc, maxs, pers))
    rownames(mos) <- numlist
    colnames(mos) <- c("total modules",
                       "singletons",
                       "% annotated modules",
                       "% characterizable features",
                       "% successful classifications",
                       "Max classes/module",
                       "90p classes/module",
                       "Max subclasses/module", 
                       "90p subclasses/module")

    logging::loginfo("Writing measures of success at candidate MMS to file.")
    write.csv(mos, file=file.path(output,"measures_of_success_by_MMS.csv"))       

    # Selecting optimum MMS and assigning features to modules
    #----------------------------------------------------------------------------------            
    # Select an optimum MMS; first MMS where pann >= 25%; fann >= 15%; hscc >= 80%
    if(is.null(min_module_size)){
      chosen.mms <- as.numeric(head(rownames(mos[which(mos[,3] >= 25 & mos[,4] >= 15 & mos[,5] >= 80),]),1))
    }else{
      chosen.mms <- as.numeric(min_module_size)
    }
    logging::loginfo(paste0("Optimum Minimum Module Size for this dataset: ", as.numeric(head(rownames(mos[which(mos[,3] >= 25 & mos[,4] >= 25 & mos[,5] >= 80),]),1)))) 
    logging::loginfo(paste0("Minimum Module Size used for this dataset (user provided): ", chosen.mms)) 
    mod.assn$module <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                               distM = as.matrix(w), 
                                                               deepSplit = TRUE, 
                                                               pamRespectsDendro = TRUE,
                                                               minClusterSize = chosen.mms,
                                                               verbose = 0))
    df.module <- as.data.frame(cbind(anno[rownames(mod.assn),], mod.assn$module))
    names(df.module)[names(df.module) == "mod.assn$module"] <- "Module"
    logging::loginfo(paste0("Total modules found: ", max(mod.assn$module)))


    logging::loginfo("Writing module assignments to file.")                                                            
    write.csv(df.module, file=file.path(output, "MACARRoN_module_assignments.csv"))


    #================================================================================== 
    # Relative abundance calculation
    #==================================================================================   
    print("Initiating relative abundance calculations")
    fint <- as.data.frame(assay(se))
    fint <- fint[rownames(mod.assn),]
    fint <- t(fint)

    # Finding anchor compound for each module
    #----------------------------------------------------------------------------------   
    anchors <- NULL
    modules <- sort(unique(mod.assn$module))
    logging::loginfo("Finding anchor metabolic features for modules") 

    find.anchor <- function(m){
      logging::logdebug(paste0("Finding anchor metabolic feature for module: ", m)) 
      if(dim(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])[1] > 0){
      #known features in the module
         r <- rownames(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])
      }else{
      #all features in the module
         r <- rownames(mod.assn[which(mod.assn$module == m),])
      }
      means <- NULL
      for (g in grps){
        ind <- se[[ptype]] == g
        if(length(r) > 1){
          gmean <- apply(fint[ind,r],2,function(x) mean(x, na.rm=TRUE))
          d <- data.frame(m, g, max(gmean), stringsAsFactors = FALSE)
          colnames(d) <- c("module","grp","mma")
        }else{
          d <- data.frame(m, g, mean(fint[ind,r], na.rm=TRUE), stringsAsFactors = FALSE)
          colnames(d) <- c("module","grp","mma")
        }
      means <- rbind(means,d)
    }
    tail(means[order(means[3]),],1)
    }
    f.list <- do.call(rbind, lapply(modules, find.anchor))
    mod.assn$feature <- rownames(mod.assn)
    anchors <- merge(mod.assn, f.list, by="module")
    rownames(anchors) <- anchors$feature
    mod.assn <- mod.assn[,-3]

    # Relative abundance
    #----------------------------------------------------------------------------------   
    logging::loginfo("Calculating relative abundance of metabolic features assigned to modules.") 
    # singletons (I)
    sing.ra <- rep(1, length(anchors[which(anchors$module==0),1]))
    sing.ra <- as.data.frame(sing.ra)
    rownames(sing.ra) <- rownames(anchors[which(anchors$module==0),])
    colnames(sing.ra) <- "rel.abun."

    # features assigned to modules (II)
    cal.relab <- function(i){
    logging::logdebug(paste0("Calculating relative abundance for feature ",i))
    getMean <- function(g){
      ind <- se[[ptype]] == g
      gmean <- mean(fint[ind,i], na.rm=TRUE)
      gmean
    }
    all.means <- sapply(grps, getMean)
    relab.i <- max(all.means)/as.numeric(as.character(anchors[i,"mma"]))
    }
    members <- rownames(anchors[which(anchors$module!=0),])
    memb.ra <- as.data.frame(do.call(rbind, lapply(members, cal.relab)))
    rownames(memb.ra) <- members
    colnames(memb.ra) <- "rel.abun."

    # combining (I) and (II)
    my.relab <- as.data.frame(rbind(sing.ra, memb.ra))
    my.relab$module <- mod.assn[row.names(my.relab),"module"]
    my.relab <- my.relab[rownames(mod.assn),]


    #================================================================================== 
    # Calculating q-value of differential abundance with MaAsLin2
    #================================================================================== 
    print("Initiating q-value calculations with MaAsLin2")

    for (lib in c('Maaslin2','plyr')) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }
    # getting input data and metadata for MaAsLin2
    fint <- as.data.frame(assay(se))
    fint <- fint[rownames(mod.assn),]
    fint[is.na(fint)] <- 0
    fint <- log2(fint + 1)
    meta <- as.data.frame(colData(se))

    # setting the control condition
    if(is.null(control_condition)){
      control_condition <- levels(as.factor(meta[,ptype]))[1]
    }else{
      control_condition <- control_condition
      meta[,ptype] <- relevel(as.factor(meta[,ptype]), ref=control_condition)
    }
    folder_name <- file.path(output, "Maaslin2_results")

    # Linear model for q-value with MaAsLin2
    #----------------------------------------------------------------------------------   
    fit.data <- Maaslin2(input_data = fint, 
                       input_metadata = meta, 
                       output = folder_name, 
                       fixed_effects = fixed_effects, 
                       random_effects = random_effects,
                       normalization = "NONE",
                       transform = "NONE",
                       min_prevalence = 0,
                       min_abundance = 0,
                       min_variance = 0,
                       cores = cores,
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)

    my.qvalue <- as.data.frame(fit.data$results[which(fit.data$results$metadata == ptype),
                              c("feature","metadata","value","coef","pval")])  
    my.qvalue <- ddply(my.qvalue, .(metadata,value), 
                              transform, qvalue=as.numeric(p.adjust(as.numeric(pval), "BH")))  

    #================================================================================== 
    # Calculating effect size of differential abundance
    #================================================================================== 
    print("Initiating effect size calculations")
    
    # calculating mean abundance in each condition
    fint <- t(fint)
    get.mean <- function(g){
      logging::loginfo(paste0("Calculating mean abundances of metabolic features in condition: ",g))
      ind <- se[[ptype]] == g
      m <- sapply(row.names(mod.assn), function(f) mean(fint[ind,f]))
    }
    all.means <- as.data.frame(do.call(cbind, lapply(grps, get.mean)))
    colnames(all.means) <- grps

    # calculating effect size
    get.es <- function(c){
      logging::loginfo(paste0("Calculating effect sizes of metabolic features in condition: ",c))
      es <- all.means[,c]- all.means[,control_condition]
    }

    # identifying non-control conditions in <metadata variable>
    test.grps <- setdiff(grps, control_condition)
    my.effsize <- as.data.frame(do.call(cbind, lapply(test.grps, get.es)))
    colnames(my.effsize) <- test.grps
    rownames(my.effsize) <- rownames(mod.assn)

    #================================================================================== 
    # Integrating ranks and prioritizing bioactives in each non-control condition
    #================================================================================== 
    for (lib in c('psych')) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }
    print("Initiating prioritization")
    all_prioritized <- NULL
    for (t in test.grps){
      # collect all scores
      logging::loginfo(paste0("Prioritizing metabolic features in condition: ",t))
      qval <- my.qvalue[which(my.qvalue$value == t),]
      rownames(qval) <- qval$feature
      qval <- qval[rownames(mod.assn),]
      all_param <- as.data.frame(cbind(my.relab[,"rel.abun."],my.effsize[,t],qval[,"qvalue"]))
      colnames(all_param) <- c("relab","efs","qval")
      rownames(all_param) <- rownames(mod.assn)

      # direction of perturbation and annotations
      all_param$perturbation <- ""
      all_param[which(all_param$efs < 0),"perturbation"] <- paste0("depleted in ", t)
      all_param[which(all_param$efs > 0),"perturbation"] <- paste0("enriched in ", t)
      all_param$efs <- abs(all_param$efs)
      all_param$module <- mod.assn[row.names(all_param),"module"]
      all_param$standard_identifier <- mod.assn[row.names(all_param),1]
      all_param$association <- 0
      all_param[which(all_param$standard_identifier != ""),"association"] <- 1

      # deprioritizing features that are not associated with a known metabolite
      annotated_modules <- unique(all_param[which(all_param$standard_identifier != "" & all_param$module != 0),"module"])
      all_param[which(all_param$module %in% annotated_modules),"association"] <- 1

      # ranks
      logging::loginfo("Assigning relative abundance ranks")
      all_param$relab_rank <- rank(all_param$relab)
      logging::loginfo("Assigning effect size ranks")
      all_param$efs_rank <- rank(all_param$efs)
      logging::loginfo("Assigning q-value ranks")
      all_param$qval_rank <- rank(-all_param$qval)
      #logging::loginfo("Assigning association with a standard metabolite ranks")
      #all_param$assoc_rank <- rank(all_param$association)

      # meta-rank
      logging::loginfo("Integrating ranks and prioritizing bioactives")
      all_param$meta_rank <- harmonic.mean(t(all_param[,8:10]))
      ranked_features <- all_param[order(-all_param$meta_rank),]
      ranked_features$priority <- sapply(seq(1:nrow(ranked_features)), function(p) p*100/(nrow(ranked_features)))
      colnames(ranked_features) <- c("relative_abundance",
                                     "effect_size",
                                     "q_value",
                                     "perturbation",
                                     "module",
                                     "standard_identifier",
                                     "associated_with_standard",
                                     "relative_abundance_rank",
                                     "effect_size_rank",
                                     "q_value_rank",
                                     "meta_rank",
                                     "priority_percentage")
      final_df <- as.data.frame(cbind(ranked_features, anno[rownames(ranked_features),]))
      file_name <- paste0("prioritized_bioactives_",t,".csv")
      file_loc = file.path(output, file_name)
      logging::loginfo(paste0("Writing prioritized bioactives in ",t," to file: ",file_loc))
      write.csv(final_df, file=file_loc)

      char_ranked <- ranked_features[which(ranked_features$associated_with_standard == 1),]
      char_df <- as.data.frame(cbind(char_ranked, anno[rownames(char_ranked),]))
      file_name <- paste0("characterizable_bioactives_",t,".csv")
      file_loc = file.path(output, file_name)
      logging::loginfo(paste0("Writing characterizable bioactives in ",t," to file: ",file_loc))
      write.csv(char_df, file=file_loc)

      #Writing Summarized Input file ()
      logging::loginfo(paste0("Writing summarized input file to ",output,"/summarized_input.csv"))
      sum.input <- as.data.frame(cbind(se[[ptype]],t(assay(se))))
      names(sum.input)[names(sum.input) == 'V1'] <- ptype
      write.csv(sum.input, file = file.path(output,"summarized_input.csv"))  
      all_prioritized <- rbind(all_prioritized, final_df)
    }
    return(all_prioritized)                                                                                                                                                                                                       
  }


if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  # check three positional arguments are provided
  if (length(positional_args) != 4) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments",
            )
    )
  }
  se <- MACARRoN(positional_args[1],
               positional_args[2],
               positional_args[3],
               positional_args[4],
               current_args$output,
               current_args$metadata_variable,
               current_args$control_condition,
               current_args$min_prevalence,
               current_args$execution_mode,
               current_args$standard_identifier,
               current_args$min_module_size,
               current_args$fixed_effects,
               current_args$random_effects,
               current_args$reference,
               current_args$cores
               )     
}
