################################################################################
##
## Importing results from contaminant single org analysis
## 
################################################################################
contamSingleOrgResults <- list.dirs("data/sim_contam",
                              full.names = TRUE, 
                              recursive = FALSE)  %>%
      grep(pattern = "uid", value = TRUE) %>% 
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0)  %>% # excluding directories without results
      map_df(read_tsv,skip = 1,.id = "Query")

cache("contamSingleOrgResults")
