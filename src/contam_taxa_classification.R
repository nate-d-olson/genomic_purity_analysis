library(ProjectTemplate)
load.project()
################################################################################
##
## Taxize classification for target, contaminant and match taxids
## 
################################################################################

get_taxid_classification_contam <- function(queryIdTbl, contamResults, taxidClassification){
      ## get taxid for query and contam orgs
      contam_uid <- with(contamResults, c(target, contaminant)) %>% unique()
      query_taxid <-  filter(queryIdTbl, query_uid %in% contam_uid) %>% 
            .$query_taxid %>% unique()
      
      match_genome <- contamResults %>%
            filter(!(Genome %in%  # filtering unknown/ unclassified taxids
                           c("ti|12908|org|unclassified_sequences",
                             "ti|0|org|Unknown."))) %>% .$Genome
      
      match_taxid <- str_match(string = match_genome,
                               pattern = "ti\\|(.*)\\|org")[,2] %>%
            unique()
      
      contam_taxids <- c(query_taxid, match_taxid) 

      ## exclude taxids already in taxidClassification
      taxids_single <- names(taxidClassification)
      search_taxid <- contam_taxids[!(contam_taxids %in% taxids_single)]
      
      if(length(search_taxid) > 0){
            taxidClassification <- search_taxid %>% 
                  classification(db = 'ncbi') %>% 
                  c(taxidClassification, .)
            
      }else{
            print("No new taxids")
      }
      taxidClassification
}


taxidClassification <-  get_taxid_classification_contam(queryIdTbl, 
                                                        contamResults, 
                                                        taxidClassification)


## Adding to project cache
cache("taxidClassification")