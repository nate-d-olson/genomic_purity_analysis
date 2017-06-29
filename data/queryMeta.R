################################################################################
##
## Extracting 
## 
################################################################################

parse_rpt <- function(file_name){
      read_lines(file_name) %>% str_replace(" =", ":") %>% 
            data_frame(rpt = .) %>% 
            separate(rpt, c("var","value"), ": ",extra = "merge") %>% 
            filter(var %in% c("Accession","GI","DNA  length","Taxname","Taxid")) %>% 
            spread(var, value)
}

parse_rpt_set <- function(file_names){
      file_names %>% map_df(parse_rpt)
}

queryMeta <- list.dirs("data/Bacteria",
                              full.names = TRUE, 
                              recursive = FALSE)  %>%
      set_names(basename(.)) %>%
      map(list.files, pattern = "*rpt",full.names = TRUE) %>%
      map_df(parse_rpt_set,.id = "Query") %>% 
      ## Accession, GI, and length not in rpt file for uid215084
      mutate(Accession = if_else(Query == "Escherichia_coli_c321D_uid215084", "CP006698.1", Accession),
             GI = if_else(Query == "Escherichia_coli_c321D_uid215084", "549811571", GI),
             `DNA  length` = if_else(Query == "Escherichia_coli_c321D_uid215084", "4643553", `DNA  length`))

 cache("queryMeta")