################################################################################
##
## Table with uid, gi, and taxid for queries
## 
################################################################################

# ## extract gi and uid from singleOrgRef
# make_id_tbl <- function(singleOrgRef){
#       single_org_id <- singleOrgRef %>%
#             mutate(uid = str_match(file, "_uid(.*)/")[,2],
#                    gi = str_match(seq_id, "gi\\|(.*)\\|[gb,emb,dbj]")[,2])
#       
#       query_id_tbl <- single_org_id %>% select(uid, gi)
#       
#       ## get taxid for gi
#       query_gi <- query_id_tbl$gi %>% unique()
#       query_taxid <- query_gi %>% map(~genbank2uid(id = . )) %>% map(1) %>% 
#             flatten_chr()
#       
#       ## joining all three
#       query_id_tbl <- data_frame(gi = query_gi, taxid = query_taxid) %>% 
#             left_join(query_id_tbl, .)
#       
#       ## get taxid for gi's with na
#       taxid_na <- query_id_tbl %>% filter(is.na(taxid))
#       query_id_na <- single_org_id %>% right_join(taxid_na)
# 
#       na_taxid_taxa <- query_id_na$file %>% dirname() %>% basename() %>%
#             str_replace("_"," ") %>% str_replace("_.*","")
# 
#       ## taxid based on species names
#       na_taxid <- data_frame(uid = query_id_na$uid,
#                         gi = query_id_na$gi,
#                              taxid = get_ids(na_taxid_taxa, db = "ncbi")$ncbi)
#       query_id_tbl %>% filter(!is.na(taxid)) %>% 
#             bind_rows(na_taxid) %>% select(-gi) %>% unique()
# }
# 
# query_id_tbl <- make_id_tbl(singleOrgRef)
