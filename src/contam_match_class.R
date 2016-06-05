library(ProjectTemplate)
load.project()
################################################################################
##
## Classification contamination matches
## 
################################################################################

contam_match_df <- contamMixResults %>% 
      filter(!(Genome %in%  # filtering unknown/ unclassified taxids
                     c("ti|12908|org|unclassified_sequences", 
                       "ti|0|org|Unknown."))) %>%
      mutate(match_taxid = str_match(string = Genome, 
                                     pattern = "ti\\|(.*)\\|org")[,2]) %>% 
      select(contam_ds, target, contaminant, mix, match_taxid, `Final Guess`) %>% 
      left_join(queryIdTbl, by = c("target"="query_uid")) %>% 
      rename(target_taxid = query_taxid, target_name = Query) %>% 
      left_join(queryIdTbl, by = c("contaminant"="query_uid")) %>% 
      rename(contam_taxid = query_taxid, contam_name = Query)     

contam_match_taxid <- contam_match_df %>% 
      select(match_taxid, contam_taxid, target_taxid) %>% distinct()

### Finding Lowest Common Match Level
contamMixMatch <- data_frame()
for(i in 1:nrow(contam_match_taxid)){
      match_taxid <- contam_match_taxid$match_taxid[i]
      match_class <- taxidClassification[[match_taxid]]
      
      target_taxid <- contam_match_taxid$target_taxid[i]
      target_class <- taxidClassification[[target_taxid]]

      contam_taxid <- contam_match_taxid$contam_taxid[i]
      contam_class <- taxidClassification[[contam_taxid]]

      target_match_class_df <- inner_join(match_class, target_class) %>% 
            filter(rank != "no rank") %>% tail(1) %>% 
            rename(lca_target_name = name, 
                   lca_target_rank = rank, 
                   lca_target_id = id) %>% 
            mutate(target_taxid = target_taxid,
                   match_taxid = match_taxid)
      
      contam_match_class_df <- inner_join(match_class, contam_class) %>% 
            filter(rank != "no rank") %>% tail(1) %>% 
            rename(lca_contam_name = name, 
                   lca_contam_rank = rank, 
                   lca_contam_id = id) %>% 
            mutate(contam_taxid = contam_taxid,
                   match_taxid = match_taxid)
      
      contamMixMatch <- bind_cols(target_match_class_df, contam_match_class_df) %>% 
            bind_rows(contamMixMatch, .)
}

## Adding to project cache
cache("contamMixMatch")

contamMixMatchResults <- contamMixMatch %>% left_join(contam_match_df)

## Check total Final guess for each dataset is greater than 1, 42 entries greater than 1 but due to rounding
# contamMixMatchResults %>% group_by(contam_ds) %>% summarise(total_guess = sum(`Final Guess`)) %>% 
#      mutate(check_round = total_guess - 1) %>% filter(check_round > 0)

## Adding to project cache
cache("contamMixMatchResults")
