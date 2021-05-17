#!/usr/bin/env Rscript

###############################################################################

# MACARRoN Utility: Assign chemical taxonomy to annotated compounds using HMDB ID or PubChem CID

###############################################################################

#load the required libraries, report an error if they are not installed

for (lib in c('data.table', 'xml2', 'RJSONIO', 'RCurl')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

args <- list()
args$input_annotations <- NULL

###############################################################################
# Add command line arguments 
###############################################################################

options <-
  optparse::OptionParser(usage = paste(
    "%prog [Inputs]",
    " <feature_annotations.csv> "
    )
)

###############################################################################
# Main decorateID function
###############################################################################

decorateID <- function(
    input_annotations
)
{
    # Read feature annotation file
    if (is.character(input_annotations)) {
      feat_anno <- read.csv(input_annotations, row.names=1)
      } else {
      feat_anno <- input_annotations
    }

    ID_list <- unique(feat_anno[,1][feat_anno[,1] != ""])

    # function if HMDB Accession
    ChemTax_HMDB <- function(h){
        hmdb_url <- paste0("https://hmdb.ca/metabolites/",h,".xml")
        Sys.sleep(0.1)
        if(url.exists(hmdb_url)){
          hmdb_page <- read_xml(hmdb_url) 
          node_class <- xml_text(xml_find_all(hmdb_page, "//class"))
          node_subclass <- xml_text(xml_find_all(hmdb_page, "//sub_class"))
          cbind(h, node_subclass, node_class)
        }
    }

    # function if PubChem Compound ID
    ChemTax_PubChem <- function(p){
        pubchem_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",p,"/property/InChIKey/TXT")
        inchi_txt <- RCurl::getURL(pubchem_url)
        inchi_key <- gsub('.{1}$','',inchi_txt)
        classy_url <- paste0("http://classyfire.wishartlab.com/entities/",inchi_key,".json")
        Sys.sleep(0.1)
        if(url.exists(classy_url)){
          classy_page <- fromJSON(classy_url)
          node_class <- as.character(classy_page[['class']][1])
          node_subclass <- as.character(classy_page[['subclass']][1])
          cbind(p, node_subclass, node_class)
        }
    }

    # create the final dataframe
    if(ID_list[1] %like% "HMDB"){
        dat <- do.call(rbind, lapply(ID_list, function(x) ChemTax_HMDB(x)))
    }else{
        dat <- do.call(rbind, lapply(ID_list, function(x) ChemTax_PubChem(x)))
    }
    colnames(dat) <- c(names(feat_anno)[1],"Sub_Class","Class")
    file_name <- paste0(names(feat_anno)[1],"_taxonomy.csv")
    write.csv(dat, file=file_name, row.names=FALSE)
    return(dat)
}
if (identical(environment(), globalenv()) &&
    !length(grep("^source\\(", sys.calls()))) {
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  # check three positional arguments are provided
  if (length(positional_args) < 1) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments",
            )
    )
  }
  taxonomy_table <- decorateID(positional_args[1])
}