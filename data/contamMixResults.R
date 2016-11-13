################################################################################
##
## Importing results from simulated contamination results
## 
################################################################################
contamMixResults <- list.dirs("data/sim_contam",
                              full.names = TRUE, 
                              recursive = TRUE)  %>%
      grep(pattern = "uid", invert = TRUE, value = TRUE) %>% 
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0)  %>% # excluding directories without results
      map_df(read_tsv,skip = 1,.id = "contam_ds") %>% 
      separate(contam_ds, c("target","contaminant"),sep = "-", remove = FALSE) %>% 
      separate(contaminant, c("contaminant","mix"), sep = "_")
#cache("contamMixResults")

contamSingleOrgResults <- list.dirs("data/sim_contam",
                                    full.names = TRUE, 
                                    recursive = FALSE)  %>%
      grep(pattern = "uid", value = TRUE) %>% 
      set_names(basename(.)) %>%
      map(list.files, pattern = "*report.tsv",full.names = TRUE) %>%
      keep(~length(.) > 0)  %>% # excluding directories without results
      map_df(read_tsv,skip = 1,.id = "Query")

#cache("contamSingleOrgResults")

## Combinding Contam Single Org and Mix
uid <- contamMixResults$target %>% unique()
contam_uid_df <- tibble(target = uid, contaminant = rep(list(uid), length(uid))) %>% 
      unnest() %>% 
      filter(target != contaminant)

contamMixResults <- contamSingleOrgResults %>% 
      group_by(Query) %>% 
      nest() %>% 
      mutate(target = str_match(Query, "_uid(.*)")[,2], 
             mix = "1.0") %>% 
      left_join(contam_uid_df) %>% 
      mutate(contam_ds = paste(target,"-",contaminant,"_", mix)) %>% 
      unnest() %>% 
      bind_rows(contamMixResults)

cache("contamMixResults")