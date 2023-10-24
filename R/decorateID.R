#' Create a chemical taxonomy table for annotated metabolic features.
#' 
#' @param input_annotations a dataframe (features x annotations) containing the available feature annotations.
#' ^^Column 1 must contain standard annotations such as HMDB ID or PubChem CID for 
#' the subset of identified/annotated metabolic features. 
#' 
#' @return tax_df input_taxonomy-dataframe containing ID (HMDB or PubChem), chemical sub class and chemical class of annotated metabolic features.
#' 
#' @examples 
#' prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
#' annotations_df = read.csv(file = prism_annotations, row.names = 1)
#' input_taxonomy <- decorateID(annotations_df)
#' 
#' @export


decorateID <- function(input_annotations)
{
    ID_list <- unique(input_annotations[,1][input_annotations[,1] != ""])
    
    # function if HMDB Accession
    ChemTax_HMDB <- function(h){
        hmdb_url <- paste0("https://hmdb.ca/metabolites/",h,".xml")
        if(length(grep("Location", curlGetHeaders(hmdb_url))) > 0){
          h_new <- stringr::str_extract(grep("Location", curlGetHeaders(hmdb_url), value = TRUE), "HMDB[0-9]*")
          hmdb_url <- paste0("https://hmdb.ca/metabolites/",h_new,".xml")
        }
        Sys.sleep(sample(5,1)*0.1)
        if(isFALSE(httr::http_error(hmdb_url))){
          hmdb_page <- xml2::read_xml(hmdb_url) 
          node_class <- xml2::xml_text(xml2::xml_find_all(hmdb_page, "//class"))
          node_subclass <- xml2::xml_text(xml2::xml_find_all(hmdb_page, "//sub_class"))
          cbind(h, node_subclass, node_class)
        }
    }

    # function if PubChem Compound ID
    ChemTax_PubChem <- function(p){
        pubchem_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",p,"/property/InChIKey/TXT")
        options(warn=-1)
        inchi_txt <- httr::content(httr::GET(pubchem_url), as="text")
        inchi_key <- gsub('.{1}$','',inchi_txt)
        classy_url <- paste0("http://classyfire.wishartlab.com/entities/",inchi_key,".json")
        Sys.sleep(sample(5,1)*0.1)
        if(isFALSE(httr::http_error(classy_url))){
          classy_page <- RJSONIO::fromJSON(classy_url)
          node_class <- as.character(classy_page[['class']][1])
          node_subclass <- as.character(classy_page[['subclass']][1])
          cbind(p, node_subclass, node_class)
        }
    }

    # create the final dataframe
    if(grepl("HMDB", ID_list[1])){
        tax_df <- do.call(rbind, lapply(ID_list, function(x) ChemTax_HMDB(x)))
    }else{
        tax_df <- do.call(rbind, lapply(ID_list, function(x) ChemTax_PubChem(x)))
    }
    colnames(tax_df) <- c(names(input_annotations)[1],"Sub_Class","Class")
    tax_df <- as.data.frame(tax_df)
    tax_df
}
