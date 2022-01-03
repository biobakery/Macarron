#' Create a chemical taxonomy table for annotated metabolic features.
#' 
#' @param feat_anno a dataframe (features x annotations) containing the available feature annotations.
#' ^^Column 1 must contain standard annotations such as HMDB ID or PubChem CID for 
#' the subset of identified/annotated metabolic features. 
#' 
#' @return dataframe containing ID (HMDB or PubChem), chemical sub class and chemical class of annotated metabolic features.
#' 
#' chem_tax <- decorateID(feat_anno)
#' 
#' @export


decorateID <- function(feat_anno)
{
    ID_list <- unique(feat_anno[,1][feat_anno[,1] != ""])
    
    for (lib in c('data.table', 'xml2', 'RJSONIO', 'RCurl')) {
      suppressPackageStartupMessages(require(lib, character.only = TRUE))
    }

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
        tax_df <- do.call(rbind, lapply(ID_list, function(x) ChemTax_HMDB(x)))
    }else{
        tax_df <- do.call(rbind, lapply(ID_list, function(x) ChemTax_PubChem(x)))
    }
    colnames(tax_df) <- c(names(feat_anno)[1],"Sub_Class","Class")
    file_name <- paste0(names(feat_anno)[1],"_taxonomy.csv")
    write.csv(tax_df, file=file_name, row.names=FALSE)
    return(tax_df)
}