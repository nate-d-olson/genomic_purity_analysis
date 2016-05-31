## Methods for obtaining taxa lineage and matches using NCBI taxa sqlite database
library('ProjectTemplate')
load.project()

## NCBI Taxa DB
ncbi_taxid_db <- "/Volumes/Transcend/genomic_purity/data/ncbi_taxa_db/ncbi_taxa.sqlite"
ncbi_tax_id <- src_sqlite(path = ncbi_taxid_db)

gi_id <- tbl(ncbi_tax_id,"gi_taxid")
ncbi_nodes <- tbl(ncbi_tax_id,"nodes")
ncbi_names <- tbl(ncbi_tax_id,"names")

## Taxid and GI for query and matches in database
# query gi
query_gi <- singleOrgRef$gi %>% unique()

# query taxid
query_gi_id <- gi_id %>% filter(gi %in% query_gi) %>% collect() %>% distinct()
query_taxid <- query_gi_id$taxid %>% unique()

# match taxid
match_taxid <- single_org_dat$match_taxid %>% unique()


## Taxonomic Table
tax_ids <- c(query_taxid, match_taxid) %>% unique()
taxonomy_table <- data_frame()
while(length(tax_ids) > 1){
      temp_tbl <- ncbi_nodes %>% 
            filter(tax_id %in% tax_ids) %>% 
            select(tax_id, parent_id, rank) %>% 
            collect()
      tax_ids <- temp_tbl$parent_id %>% unique()
      taxonomy_table <-  bind_rows(taxonomy_table, temp_tbl) %>% 
            unique()
}

## Lineage Dataframe
get_lineage <- function(query_id, tax_tbl = taxonomy_table){
      query_lineage <- list()
      not_root <- TRUE
      query_lvl <- 0
      current_id <- query_id
      query_lin <- data_frame()
      
      if(!(query_id %in% taxonomy_table$tax_id)){
            print(paste("Tax ID: ", query_id," is not in taxonomy_table"))
            return(query_lin)
      }
      while(not_root){
            current_lvl <- taxonomy_table %>% filter(tax_id == current_id) %>% 
                  mutate(query_id = query_id, query_lvl = query_lvl)
            
            query_lin <- query_lin %>%  
                  bind_rows(current_lvl)
            
            if(current_id == current_lvl$parent_id){
                  not_root <- FALSE
            }else{
                  current_id <- current_lvl$parent_id
                  query_lvl <- query_lvl + 1
            }
      }
      query_lin
}

tax_ids <- c(query_taxid, match_taxid) %>% unique()
tax_lineage <- tax_ids %>% map_df(get_lineage)

tax_lineage <-  query_gi_id %>% mutate(query_id = as.character(taxid)) %>% 
      right_join(tax_lineage)


## The following tax_ids are not in the ncbi db
## data_frame(tax_id = tax_ids) %>% anti_join(taxonomy_table)
##  tax_id
## (chr)
## 1      0
## 2 145481
## 3 544436
## 4 508765