################################################################################
##
## query id table with query name, uid, and taxid
## 
################################################################################

make_query_id_tbl <- function(singleOrgResults, queryMeta){
      results_query <- singleOrgResults$Query %>% unique()
      queryMeta %>% filter(Query %in% results_query) %>% 
            select(Query, Taxid) %>% unique() %>% 
            ## removing duplicate Taxid - see queryMeta_check.Rmd
            filter(!(Query == "Clostridium_perfringens_SM101_uid12521" & 
                           Taxid == "396359"),
                   !(Query == "Listeria_monocytogenes_serotype_1_2b_SLCC2755_uid50379" &
                           Taxid == "1639")) %>% 
            mutate(query_uid = str_match(Query, "_uid(.*)")[,2]) %>% 
            rename(query_taxid = Taxid)
            
}

queryIdTbl <- make_query_id_tbl(singleOrgResults, queryMeta)

cache("queryIdTbl")