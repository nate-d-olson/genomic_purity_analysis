## Lowest common ancestor using ncbi taxa db sqlit source
uniq_seq_gi <- unique(singleOrgRef$gi)

gi_tbl <- gi_id %>% filter(gi %in% uniq_seq_gi) %>% collect()

uniq_seq_id <- singleOrgRef %>% left_join(gi_tbl)

## Compare taxa
ref_uid2taxid <- uniq_seq_id %>% select(uid, taxid)

# removing duplicate entries for multiple gi
tax_lineage_uniq <- tax_lineage %>% select(-gi) %>% distinct()

match_results <- single_org_dat %>% select(query_uid, match_taxid) %>% 
      mutate(query_uid = as.integer(query_uid),
             match_taxid = as.integer(match_taxid)) %>% 
      left_join(ref_uid2taxid,by = c("query_uid"="uid")) %>% 
      rename(query_taxid = taxid)

# get_match_level <- function(query_taxid, match_taxid){
#       if(!(query_taxid %in% lin_df$taxid)){
#             return("missing_query_id")
#       }
#       if(!(match_taxid %in% lin_df$taxid)){
#             return("missing_match_id")
#       }
#       query_lin <- lin_df %>% 
#             filter(taxid == query_taxid)
#       match_lin <- lin_df %>% 
#             filter(taxid == match_taxid)
#       
#       for(lvl in c("species","genus","family","order","class","phylum")){
#             if(!is.na(query_lin[[lvl]]) & !is.na(match_lin[[lvl]])){
#                   if(query_lin[[lvl]] == match_lin[[lvl]]){
#                         return(lvl)
#                   }   
#             }
#             
#       }
#       "no match"
# }
get_match_level <- function(query_taxid, match_taxid){
      query_lineage <- tax_lineage_uniq %>% 
            filter(query_id == query_taxid) %>% 
            select(-taxid) 
      match_lineage <- tax_lineage_uniq %>% 
            filter(query_id == match_taxid) %>% 
            select(-taxid) %>% rename(match_id = query_id,
                                      match_lvl = query_lvl) 
      semi_join(query_lineage, match_lineage) %>% 
            top_n(1, wt = -query_lvl) %>% .$rank
}


# Missing query_taxid
match_results %>% filter(is.na(query_taxid)) %>% .$query_uid %>% unique()

uniq_match <- match_results %>% select(query_taxid, match_taxid) %>% 
      distinct()
match_lvl <- c()
for(i in 1:nrow(uniq_match)[1:20]){
      query_taxid <- uniq_match$query_taxid[i]
      match_taxid <- uniq_match$match_taxid[i]
      if(is.na(query_taxid)){
            lvl <- "q_na"
      }else if(query_taxid == match_taxid){
            lvl <- "exact"
      }else{
            lvl <- get_match_level(query_taxid, match_taxid)
            if(length(lvl) == 0){
                  lvl <- "no match"
            }      
      }
      match_lvl <- c(match_lvl, lvl)
}

uniq_match$match_lvl <- match_lvl

pathoid_match <- single_org_dat %>% 
      mutate(match_taxid = as.integer(match_taxid),
             query_uid = as.integer(query_uid)) %>% 
      left_join(uniq_match)


      
## Issue with multiple entries????
# pathoid_match %>% distinct() %>% group_by(org, match_lvl) %>% 
#       filter(!(match_lvl %in% c("missing_match_id","missing_query_id"))) %>% 
#       summarise(total_guess = sum(`Final Guess`)) %>% filter(total_guess > 1)
# 
# pathoid_match %>% filter(org == "Bacillus_thuringiensis_serovar_IS5056_uid187142", match_lvl == "species") %>% View()


## Most matches have missing query or match taxonomic ids
ggplot(pathoid_match) + geom_bar(aes(x = match_lvl))

pathoid_clean <- pathoid_match %>%
      select(org, match_taxid, match_lvl, `Final Guess`) %>% 
      group_by(org, match_lvl) %>% distinct() %>% 
      summarise(total_guess = sum(`Final Guess`))

pathoid_clean %>% 
      ggplot() + 
            geom_boxplot(aes(x = match_lvl, y = total_guess)) +
            theme_bw()


