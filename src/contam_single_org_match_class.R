library(ProjectTemplate)
load.project()
################################################################################
##
## Classification single org contaminant results
## 
################################################################################


contam_single_org_match <- contamSingleOrgResults %>% 
      filter(!(Genome %in%  # filtering unknown/ unclassified taxids
                     c("ti|12908|org|unclassified_sequences", 
                       "ti|0|org|Unknown."))) %>%
      mutate(match_taxid = str_match(string = Genome, 
                                     pattern = "ti\\|(.*)\\|org")[,2]) %>% #,
      # query_uid = str_match(Query, "_uid(.*)")[,2]) %>% 
      select(match_taxid, Query,`Final Guess`) %>% left_join(queryIdTbl) 

contam_single_org_match_taxid <- contam_single_org_match %>% 
      select(match_taxid, query_uid, query_taxid)

### Finding Lowest Common Match Level
match_class_df <- data_frame()
for(i in 1:nrow(contam_single_org_match_taxid)){
      match_taxid <- contam_single_org_match_taxid$match_taxid[i]
      match_class <- taxidClassification[[match_taxid]]
      
      query_taxid <- contam_single_org_match_taxid$query_taxid[i]
      query_class <- taxidClassification[[query_taxid]]
      query_genus <- query_class$name[query_class$rank == "genus"]
      if(length(query_genus) == 0){
            query_genus = "no_genus"
      }
      match_class_df <- inner_join(match_class, query_class) %>% 
            filter(rank != "no rank") %>% tail(1) %>% 
            mutate(query_taxid = query_taxid,
                   query_genus = query_genus, 
                   match_taxid = match_taxid) %>% 
            bind_rows(match_class_df, .)
}

## Adding unique duplicates caused inflated total counts
contamSingleOrgMatch <- match_class_df %>% 
      rename(lca_name = name, lca_rank = rank, lca_id = id) %>% unique()

## Adding to project cache
cache("contamSingleOrgMatch")

contamSingleOrgMatchResults <- contamSingleOrgMatch %>% left_join(contam_single_org_match)

## Adding to project cache
cache("contamSingleOrgMatchResults")