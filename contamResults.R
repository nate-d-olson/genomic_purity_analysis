################################################################################
##
## Importing results from simulated contamination results
## 
################################################################################
contamResults <- list.dirs("data/sim_contam",
                              full.names = TRUE, 
                              recursive = TRUE)  %>%
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0)  %>% # excluding directories without results
      map_df(read_tsv,skip = 1,.id = "contam_ds") %>% 
      separate(contam_ds, c("target","contaminant"),sep = "-", remove = FALSE) %>% 
      separate(contaminant, c("contaminant","mix"), sep = "_")
cache("contamResults")

contamSingleOrgResults <- contamResults %>% filter(contam_ds == target)
cache("contamSingleOrgResults")

contamMixResults <- contamResults %>% filter(contam_ds != target)
cache("contamMixResults")
