library(ProjectTemplate)
load.project()
################################################################################
##
## Taxize classification for query and match taxids
## 
################################################################################

get_taxid_classification <- function(queryIdTbl, singleOrgResults){
      query_taxid <- queryIdTbl$query_taxid %>% unique()
      
      match_genome <- singleOrgResults %>%
            filter(!(Genome %in%  # filtering unknown/ unclassified taxids
                           c("ti|12908|org|unclassified_sequences",
                             "ti|0|org|Unknown."))) %>% .$Genome
      
      match_taxid <- str_match(string = match_genome,
                               pattern = "ti\\|(.*)\\|org")[,2] %>%
            unique()
      
      taxids <- c(query_taxid, match_taxid) %>% unique()
      
      for(i in 0:12){
            print(i)
            fst <- i * 100 + 1; lst <- i * 100 + 100
            tmp_class <- taxids[fst:lst] %>% classification(db = 'ncbi')
            taxid_classification <- c(taxid_classification,tmp_class)
            Sys.sleep(30) ## to prevent timeout error
      }
      c(taxid_classification, taxids[1301:1354] %>% classification(db = 'ncbi'))
}


taxidClassification <-  get_taxid_classification(queryIdTbl, singleOrgResults)

## Adding to project cache
cache("taxidClassification")
