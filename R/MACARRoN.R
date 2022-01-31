#' MACARRoN
#' 
#' @param input_abundances a comma-delimited file or dataframe (features x samples) containing metabolic feature intensities (abundances).
#' @param input_annotations a comma-delimited file or dataframe (features x annotations) containing available feature annotations.
#' @param input_metadata a comma-delimited file or dataframe (samples x metadata) containing sample metadata.
#' @param input_taxonomy a comma-delimited file or dataframe containing the chemical class and subclass information of annotated features.
#' @param output name of the folder where MACARRoN output files will be written. Default: "MACARRoN_output".
#' @param metadata_variable the header of the column in the metadata file that identifies the main phenotypes/conditions in the study. Default: Column 2 of metadata file.
#' @param min_prevalence prevalence threshold (percentage). Default = 0.7.
#' @param execution_mode BiocParallel execution mode. Options: "serial" or "multi" Default = "serial".
#' @param standard_identifier HMDB ID or PubChem CID. Default: Column 2 of annotation file/Column 1 of annotation dataframe.
#' @param anchor_annotation Metabolite name. Default: Column 3 of annotation file/Column 2 of annotation dataframe.
#' @param min_module_size Integer that defines the size of the smallest covariance module. Default: Cube root of number of prevalent metabolic features. 
#' @param fixed_effects Covariates for linear modeling with MaAsLin2. Default: All columns of metadata dataframe.
#' @param random_effects Random effects for linear modeling with MaAsLin2. Default: NULL.
#' @param reference Reference category (factor) in categorical metadata covariates containing three or more levels. Must be provided as a string of 'covariate,reference' semi-colon delimited for multiple covariates.
#' @param cores MaAsLin2 parameter-The number of R processes to be run in parallel.
#' 
#' @return data.frame containing metabolic features listed according to their priority (potential bioactivity) in a phenotype of interest.
#' 
#' @import SummarizedExperiment
#' @import logging
#' @import data.table
#' @import WGCNA
#' @import DelayedArray
#' @import BiocParallel
#' @import ff
#' @importFrom dynamicTreeCut cutreeDynamic
#' @import Maaslin2
#' @importFrom plyr ddply
#' @importFrom stats p.adjust
#' @importFrom psych harmonic.mean
#' 
#' 
#' @export 


MACARRoN <- 
  function(
    input_abundances, 
    input_annotations, 
    input_metadata, 
    input_taxonomy,
    output = "MACARRoN_output",
    metadata_variable = NULL,
    min_prevalence = 0.7,
    execution_mode = "serial",
    standard_identifier = NULL,
    anchor_annotation = NULL,
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
    
    
    
    # Create Summarized Experiment object
    #----------------------------------------------------------------------------------
    # Create metabolic feature IDs: F[#]
    ids <- lapply(seq_len(as.numeric(nrow(feat_anno))), function(i) {paste0("F",i)})
    # Assign newly created feature IDs as rownames for intensity and annotation tables.
    rownames(feat_int) <- ids
    rownames(feat_anno) <- ids
    
    # Match sample IDs between metadata (exp_meta) and abundance data (feat_int).
    # Removes samples that do not have features/metadata.
    rownames(exp_meta) <- exp_meta[,1]
    exp_meta <- exp_meta[,-1]
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
    for (lib in c('WGCNA', 'DelayedArray', 'BiocParallel', 'ff')) {
      requireNamespace(lib, quietly = TRUE)
    }
    
    # Abundance matrix
    mat <- DelayedArray::DelayedArray(SummarizedExperiment::assay(se))
    
    # Phenotype i.e. groups/conditions
    if(is.null(metadata_variable)){
      ptype <- names(SummarizedExperiment::colData(se))[1]
      logging::loginfo("Metadata chosen for prevalence filtering: %s", ptype)
    }else{
      ptype <- metadata_variable
      logging::loginfo("Metadata chosen for prevalence filtering: %s", ptype)
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
    # Assigning prevalent metabolic features to covariance modules
    #==================================================================================
    message("Initiating module assignments")
    
    # Create bicor-based distance matrix
    #----------------------------------------------------------------------------------
    
    # Compute correlation and distance matrices
    if(nrow(mat) <= 25000){
      .getCorMat <- function(g)
      {
        inx <- se[[ptype]] == g
        gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=min_prevalence),inx]
        tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
        rownames(tmp) <- rownames(mat)
        colnames(tmp) <- colnames(gmat)
        tmp <- DelayedArray::DelayedArray(tmp)
        tmp[rownames(gmat), colnames(gmat)] <- gmat
        tmp <- log2(tmp)
        options(warn=-1)
        cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs", quick=0.05)
        options(warn=0)
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
      inx <- se[[ptype]] == g
      gmat <- mat[which(rowMeans(!is.na(mat[,inx]))>=min_prevalence),inx]
      tmp <- matrix(as.numeric(), nrow=nrow(mat),ncol=ncol(gmat))
      rownames(tmp) <- rownames(mat)
      colnames(tmp) <- colnames(gmat)
      tmp <- DelayedArray::DelayedArray(tmp)
      tmp[rownames(gmat), colnames(gmat)] <- gmat
      tmp <- log2(tmp)
      options(warn=-1)
      cmat <- WGCNA::bicor(t(tmp), use = "pairwise.complete.obs", quick=0.05)
      options(warn=0)
      cmat[is.na(cmat)] <- 0
      mat[cmat < 0] <- 0
      message(paste0(g," cmat created"))
      assign(paste0(g,"_ff"),ff::as.ff(cmat))
      message(paste0(g,"_ff created"))
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
    requireNamespace('dynamicTreeCut', quietly = TRUE)
    
    # Construct tree 
    tree <- hclust(as.dist(w), method="average")
    logging::loginfo("Tree constructed.")
    
    # Setting the minimum module size
    if(is.null(min_module_size)){
      mms = round((nrow(w))^(1/3))
    }else{
      mms = as.numeric(as.character(min_module_size))
    }
    logging::loginfo("Minimum module size used for this dataset: %d", mms)
    
    
    anno <- as.data.frame(SummarizedExperiment::rowData(se))
    if(is.null(standard_identifier)){
      mod.assn <- as.data.frame(anno[colnames(w),1])
    }else{
      mod.assn <- as.data.frame(anno[colnames(w),standard_identifier])
    }
    rownames(mod.assn) <- colnames(w)
    
    
    # Evaluating Measures of Success
    
    # Range of mms for measures of success
    mms.list <- c(mms - 10, mms - 5, mms, mms + 5, mms + 10)
    mms.list <- mms.list[which(mms.list > 1)]
    
    # Measures of success
    sing <- NULL # singletons
    totc <- NULL # total modules
    pann <- NULL # % annotated modules
    hscc <- NULL # % features with the same annotation that are assigned to the same module
    maxc <- NULL # max classes per module
    maxs <- NULL # max subclasses per module
    perc <- NULL # 90th percentile classes per module
    pers <- NULL # 90th percentile subclasses per module
    fham <- NULL # % features in homogeneously annotated modules
    
    for (n in mms.list){
      mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                              distM = as.matrix(w), 
                                                              deepSplit = TRUE, 
                                                              pamRespectsDendro = TRUE,
                                                              minClusterSize = n))
      # singletons
      sing <- rbind(sing, length(which(mod.assn[,2] == 0)))
      
      # total modules
      totc <- rbind(totc, max(mod.assn[,2]))
      
      # % annotated modules
      ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
      ann.mod <- ann.mod[which(ann.mod > 0)]
      pann <- rbind(pann, round((length(ann.mod)*100)/max(mod.assn[,2]),2))
      
      # % features with the same annotation that are assigned to the same module
      same.mod <- NULL
      diff.mod <- NULL
      times.obs <- as.data.frame(table(mod.assn[,1][which(mod.assn[,1] != "")]))
      times.obs <- times.obs[which(times.obs$Freq > 1),]
      for (t in times.obs$Var1){
        if(length(unique(mod.assn[which(mod.assn[,1] == t),2])) > 1){
          diff.mod <- rbind(diff.mod,t)}else{
            same.mod <- rbind(same.mod,t)
          }
      }
      hscc <- rbind(hscc,round(length(same.mod)*100/(length(same.mod)+length(diff.mod)),2))
      
      # chemical homogeneity of modules
      cls <- NULL
      scl <- NULL
      mod.assn.ann <- mod.assn[which(mod.assn[,2] %in% ann.mod & mod.assn[,1] != ""),]
      for (i in unique(mod.assn.ann[,2])){
        dat <- mod.assn.ann[which(mod.assn.ann[,2] == i & mod.assn.ann[,1] != ""),]
        mod.tax <- NULL
        for (d in unique(dat[,1])){
          feat.tax <- chem_tax[which(chem_tax[,1] == d),]
          mod.tax <- rbind(mod.tax,feat.tax)
        }
        cls <- rbind(scl, length(unique(mod.tax[,3])))
        scl <- rbind(scl, length(unique(mod.tax[,2])))
      }
      maxs <- rbind(maxs, max(scl))
      maxc <- rbind(maxc, max(cls))
      pers <- rbind(pers, quantile(scl,0.9))
      perc <- rbind(perc, quantile(cls,0.9))
      
      # % features in homogeneously annotated modules
      dat <- as.data.frame(unique(mod.assn.ann))
      dat$class <- sapply(as.character(dat[,1]), 
                          function(x) as.character(chem_tax[which(chem_tax[,1] == x),3]))
      rownames(dat) <- NULL
      dat$class <- as.character(dat$class)
      dat$class[which(dat$class == "character(0)")] <- ""
      dat <- dat[which(dat$class != ""),]
      mods.with.tax <- as.data.frame(sort(unique(dat[,2])))
      names(mods.with.tax) <- "module"
      
      # function to estimate homogeneity of each module
      findClassHomo <- function(m){
        df <- dat[which(dat[,2] == m),]
        annotated.features <- nrow(df)
        class.counts <- as.data.frame(table(df$class))
        max(class.counts$Freq)/annotated.features
      }
      
      mods.with.tax$homo <- sapply(mods.with.tax$module, function(m) findClassHomo(m))
      mods.with.tax$homo <- as.numeric(as.character(mods.with.tax$homo))
      homogeneous.mods <- unique(mods.with.tax[which(mods.with.tax$homo >= 0.75),"module"])
      
      perc.feats <- round((nrow(mod.assn[which(mod.assn[,2] %in% homogeneous.mods),])*100)/nrow(mod.assn),2)
      fham <- rbind(fham, perc.feats)
    }  
    
    # Writing results to file
    mac.mos <- data.frame(cbind(mms.list, totc, sing, pann, hscc, maxc, perc, maxs, pers, fham))
    rownames(mac.mos) <- NULL
    colnames(mac.mos) <- c("Minimum module size (MMS)",
                           "Total modules",
                           "Singletons",
                           "% Annotated modules",
                           "% Consistent assignments",
                           "Max classes per module",
                           "90p classes per module",
                           "Max subclasses per module",
                           "90p subclasses per module",
                           "% Features in HAM")
    
    file_name <- "MAC_modules_Measures_of_success.csv"
    file_loc = file.path(output, file_name)
    logging::loginfo(paste0("Writing measures of success to file: ",file_loc))
    write.csv(mac.mos, file=file_loc, row.names=FALSE)
    
    # Module assignments
    mod.assn[,2] <- as.vector(dynamicTreeCut::cutreeDynamic(dendro = tree, 
                                                            distM = as.matrix(w), 
                                                            deepSplit = TRUE, 
                                                            pamRespectsDendro = TRUE,
                                                            minClusterSize = mms))
    
    if(is.null(standard_identifier)){
      colnames(mod.assn) <- c(names(anno)[1],"module")
    }else{
      colnames(mod.assn) <- c(standard_identifier,"module")
    }
    
    logging::loginfo("Total number of modules: %d", max(mod.assn[,2]))
    
    ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
    ann.mod <- ann.mod[which(ann.mod > 0)]
    ann.mod <- as.data.frame(ann.mod)
    colnames(ann.mod) <- "module"
    
    # Annotate modules with chemical composition information based on annotated features
    assignChemTax <- function(m){
      if(m > 0){
        annotated.features <- unique(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),1])
        classes <- toString(unique(chem_tax[which(chem_tax[,1] %in% annotated.features),3]))
      }else{
        classes <- ""
      }
    }
    ann.mod$classes <- as.character(sapply(ann.mod$module, function(m) assignChemTax(m)))
    ann.mod[ann.mod == "character(0)"] <- ""
    mod.assn$classes <- as.character(sapply(mod.assn$module, function(m) ann.mod[which(ann.mod$module == m),2]))
    mod.assn[mod.assn == "character(0)"] <- ""
    mod.assn$classes <- gsub(",",";",mod.assn$classes)
    mod.assn <- as.data.frame(mod.assn)
    
    
    
    #================================================================================== 
    # Abundance versus anchor (AVA) calculation
    #==================================================================================   
    message("Initiating AVA calculations")
    fint <- as.data.frame(SummarizedExperiment::assay(se))
    fint <- fint[rownames(mod.assn),]
    fint <- t(fint)
    
    # Setting the metadata
    if(is.null(ptype)){
      ptype <- names(SummarizedExperiment::colData(se))[1]
      logging::loginfo("Metadata chosen for AVA calculation: %s", ptype)
    }else{
      ptype = ptype
      logging::loginfo("Metadata chosen for AVA calculation: %s", ptype)
    }
    grps <- unique(se[[ptype]])
    
    # Identification of the anchor in each module
    anchors <- NULL
    modules <- sort(unique(mod.assn$module))
    
    # Function for anchor and the reference phenotype
    find.anchor <- function(m){
      message(paste0("Finding anchor feature for module: ",m))
      if(dim(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])[1] > 0){
        # known features in the module
        r <- rownames(mod.assn[which(mod.assn$module == m & mod.assn[,1] != ""),])
      }else{
        # all features in the module
        r <- rownames(mod.assn[which(mod.assn$module == m),])
      }
      
      # Mean abundances in phenotypes
      means <- NULL
      for (g in grps){
        ind <- se[[ptype]] == g
        if(length(r) > 1){
          gmean <- apply(fint[ind,r],2,function(x) mean(x, na.rm = TRUE))
          d <- data.frame(m, g, max(gmean), stringsAsFactors = FALSE)
          colnames(d) <- c("module","grp","mma")
        }else{
          d <- data.frame(m, g, mean(fint[ind,r], na.rm = TRUE), stringsAsFactors = FALSE)
          colnames(d) <- c("module","grp","mma")
        }
        means <- rbind(means,d)
      }
      tail(means[order(means[3]),],1)
    }
    f.list <- do.call(rbind, lapply(modules, find.anchor))
    mod.assn$feature <- rownames(mod.assn)
    anchors <- as.data.frame(merge(mod.assn, f.list, by = "module"))
    rownames(anchors) <- anchors$feature
    
    # Abundance versus anchor
    # singletons
    sing.ava <- rep(1, length(anchors[which(anchors$module == 0),1]))
    sing.ava <- as.data.frame(sing.ava)
    rownames(sing.ava) <- rownames(anchors[which(anchors$module == 0),])
    colnames(sing.ava) <- "ava"
    
    
    # module members
    cal.ava <- function(i){
      message(paste0("Calculating AVA for feature: ",i))
      getMean <- function(g){
        ind <- se[[ptype]] == g
        gmean <- mean(fint[ind,i], na.rm = TRUE)
        gmean
      }
      all.means <- sapply(grps, getMean)
      ava.i <- max(all.means)/as.numeric(as.character(anchors[i,"mma"]))
    }
    members <- rownames(anchors[which(anchors$module != 0),])
    memb.ava <- as.data.frame(do.call(rbind, lapply(members, cal.ava)))
    rownames(memb.ava) <- members
    colnames(memb.ava) <- "ava"
    
    # Combining
    ava <- as.data.frame(rbind(sing.ava, memb.ava))
    ava$module <- mod.assn[row.names(ava),"module"]
    ava <- ava[rownames(mod.assn),]
    ava <- ava[,c(2,1)]
    ava$ava <- round(ava$ava, 6)
    
    # Assign anchors
    anno <- as.data.frame(SummarizedExperiment::rowData(se))
    if(is.null(anchor_annotation)){
      anchor.anno <- colnames(anno)[2]
    }else{
      anchor.anno <- anchor_annotation
    }
    
    ann.mod <- unique(mod.assn[which(mod.assn[,1] != ""),2])
    ann.mod <- ann.mod[which(ann.mod > 0)]
    ann.mod <- as.data.frame(ann.mod)
    colnames(ann.mod) <- "module"
    
    assignAnchor <- function(m){
      anchor.feature <- rownames(ava[which(ava$module == m & ava$ava == 1),])
      anchor.name <- as.character(anno[anchor.feature, anchor.anno])
    }
    ann.mod$anchor <- as.character(sapply(ann.mod$module, function(m) assignAnchor(m)))
    ann.mod[ann.mod == "character(0)"] <- ""
    ava$anchor <- as.character(sapply(ava$module, function(m) ann.mod[which(ann.mod$module == m),2]))
    ava[ava == "character(0)"] <- ""
    mac.ava <- as.data.frame(ava)
    
    
    #================================================================================== 
    # Calculating q-value of differential abundance with MaAsLin2
    #================================================================================== 
    message("Initiating q-value calculations with MaAsLin2")
    
    for (lib in c('Maaslin2','plyr','stats')) {
      requireNamespace(lib, quietly = TRUE)
    }
    # getting input data and metadata for MaAsLin2
    fint <- as.data.frame(SummarizedExperiment::assay(se))
    fint <- fint[rownames(mod.assn),]
    fint[is.na(fint)] <- 0
    fint <- log2(fint + 1)
    meta <- as.data.frame(SummarizedExperiment::colData(se))
    
    folder_name <- file.path(output, "Maaslin2_results")
    logging::logdebug(paste0("Writing MaAsLin2 results to folder: ",folder_name))
    
    # Linear model for q-value with MaAsLin2
    #----------------------------------------------------------------------------------   
    fit.data <- Maaslin2::Maaslin2(input_data = fint, 
                         input_metadata = meta, 
                         output = folder_name, 
                         fixed_effects = fixed_effects, 
                         random_effects = random_effects,
                         reference = reference,
                         normalization = "NONE",
                         transform = "NONE",
                         min_prevalence = 0,
                         min_abundance = 0,
                         min_variance = 0,
                         cores = cores,
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE)
    
    results_dat <- as.data.frame(fit.data$results[which(fit.data$results$metadata == ptype),
                                                  c("feature","metadata","value","coef","pval")]) 
    adjusted_dat <- plyr::ddply(results_dat, .(metadata,value), 
                          transform, qvalue=as.numeric(stats::p.adjust(as.numeric(pval), "BH")))
    rm(results_dat)                          
    mac.qval <- adjusted_dat[,c("feature","metadata","value","qvalue")]
    mac.qval                            
    
    
    #================================================================================== 
    # Calculating effect size of differential abundance
    #================================================================================== 
    
    message("Initiating effect size calculations")
    
    # calculating mean abundance in each condition
    # Abundance matrix
    fint <- as.data.frame(SummarizedExperiment::assay(se))
    fint <- fint[unique(mac.qval$feature),]
    fint[is.na(fint)] <- 0
    fint <- log2(fint + 1)
    fint <- t(fint)
    
    # Get phenotype from mac.qval
    ptype <- unique(mac.qval$metadata)
    grps <- unique(se[[ptype]])
    
    # Mean abundance of feature in each group of ptype
    get.mean <- function(g){
      logging::logdebug(paste0("Calculating mean abundances in phenotype: ",g))
      ind <- se[[ptype]] == g
      m <- sapply(colnames(fint), function(f) mean(fint[ind,f]))
    }
    all.means <- as.data.frame(do.call(cbind, lapply(grps, get.mean)))
    colnames(all.means) <- grps
    
    # Effect size
    test.grps <- unique(mac.qval$value)
    ref.grp <- setdiff(grps, test.grps)
    get.es <- function(tg){
      logging::logdebug(paste0("Calculating effect sizes in phenotype: ",tg))
      es <- all.means[,tg]- all.means[,ref.grp]
    }
    mac.es <- as.data.frame(do.call(cbind, lapply(test.grps, get.es)))
    colnames(mac.es) <- test.grps
    rownames(mac.es) <- colnames(fint)
    mac.es
    
    #================================================================================== 
    # Integrating ranks and prioritizing bioactives in each non-control condition
    #================================================================================== 
    requireNamespace('psych', quietly = TRUE)
    message("Initiating prioritization")
    
    # Case conditions
    case.grps <- unique(mac.qval$value)
    
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
      all.params$meta_rank <- psych::harmonic.mean(t(all.params[,6:8]))
      ranked.features <- all.params[order(-all.params$meta_rank),]
      rank.perc <- ecdf(all.params$meta_rank)
      ranked.features$rank_percentile <- sapply(ranked.features$meta_rank, function(x) rank.perc(x))
      ranked.features <- as.data.frame(ranked.features)
    }
    prioritized.features <- as.data.frame(do.call(rbind, lapply(case.grps, prioritize.each)))
    prioritized.features$module <- mod.assn[prioritized.features$feature,"module"]
    prioritized.features$anchor <- mac.ava[prioritized.features$feature,"anchor"]
    prioritized.features$module_composition <- mod.assn[prioritized.features$feature,"classes"]
    prioritized.features$characterizable <- 0
    prioritized.features[which(prioritized.features$anchor != ""),"characterizable"] <- "1"
    anno <- as.data.frame(SummarizedExperiment::rowData(se))
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
    
    file_name <- "prioritized_metabolites_all.csv"
    file_loc = file.path(output, file_name)
    logging::loginfo(paste0("Writing all prioritized metabolites to file: ",file_loc))
    write.csv(mac.result, file=file_loc, row.names=FALSE)
    
    file_name <- "prioritized_metabolites_characterizable.csv"
    file_loc = file.path(output, file_name)
    logging::loginfo(paste0("Writing characterizable prioritized metabolites to file: ",file_loc))
    write.csv(mac.result[which(mac.result[,9] == 1),], file=file_loc, row.names=FALSE)
    mac.result                                                                                                                                                                                                
  }
