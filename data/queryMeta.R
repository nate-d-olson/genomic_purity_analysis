################################################################################
##
## Extracting 
## 
################################################################################

parse_rpt <- function(file_name){
      read_lines(file_name) %>% str_replace(" =", ":") %>% 
            data_frame(rpt = .) %>% 
            separate(rpt, c("var","value"), ": ") %>% 
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
      map_df(parse_rpt_set,.id = "Query")


 cache("queryMeta")