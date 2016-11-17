################################################################################
##
## Importing mapping summary from contaminant single org analysis
## 
################################################################################
contamSingleOrgMapResults <- list.dirs("data/sim_contam",
                                       full.names = TRUE, 
                                       recursive = FALSE)  %>%
      grep(pattern = "uid", value = TRUE) %>% 
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0) %>% # excluding directories without results
      map_df(read_tsv,
             col_names = c("H1","aligned_reads","H2","mapped_genomes"),
             n_max = 1,.id = "Query") %>% 
      select(-H1, -H2)

cache("contamSingleOrgMapResults")