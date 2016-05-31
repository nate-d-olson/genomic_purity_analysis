################################################################################
##
## Importing results from single org analysis
## 
################################################################################
singleOrgResults <- list.dirs("data/simdir",
                              full.names = TRUE, 
                              recursive = FALSE)  %>%
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0)  %>% # excluding directories without results
      map_df(read_tsv,skip = 1,.id = "Query")

cache("singleOrgResults")