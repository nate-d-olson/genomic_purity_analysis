## Genomic Purity Publication Variables, tables, and plots
library(xtable)
library(ProjectTemplate)

# load project
wd <- getwd()
setwd("../")
load.project()
setwd(wd)

##
## single genome data ---------------------------------------------------------
##

single_tbl_df <- singleOrgMatchResults %>%
      filter(query_genus %in% c("Bacillus", "Clostridium", "Escherichia",
                                "Francisella", "Listeria", "Pseudomonas",
                                "Salmonella", "Shigella", "Staphylococcus",
                                "Yersinia")) %>%
      select(query_genus, Query, query_taxid) %>%
      rename(Taxid = query_taxid) %>% left_join(queryMeta) %>%
      distinct() %>%
      group_by(query_genus, Query) %>%
      mutate(`DNA  length` = as.numeric(`DNA  length`)) %>%
      summarise(genome_size = sum(`DNA  length`)) %>%
      group_by(query_genus) %>%
      summarise(count = n(), size_median = median(genome_size, na.rm = T)/1e6,
                size_min = min(genome_size, na.rm = T)/1e6,
                size_max = max(genome_size, na.rm = T)/1e6) %>%
      mutate(`Genome Size (Mb)` = paste0(round(size_median,2) ," (",
                                         round(size_min,2), "-",
                                         round(size_max, 2), ")"),
             query_genus = paste0("\\textit{",query_genus,"}")) %>%
      select(-size_min, -size_max, -size_median) %>%
      rename(Genus = query_genus, N = count) %>% arrange(desc(N))


species_lvls <- c("species group", "species subgroup", "subspecies")
single_org_cum <- singleOrgMatchResults %>%
      filter(query_genus %in% c("Bacillus","Clostridium","Escherichia",
                                "Francisela","Listeria","Pseudomonas",
                                "Salmonella","Shigella","Staphylococcus","Yersinia")) %>%
      mutate(lca_rank = if_else(lca_rank %in% species_lvls,
                                "species", lca_rank)) %>%
      group_by(Query, lca_rank, query_genus) %>%
      summarise(total_prop = sum(`Final Guess`)) %>%
      spread(lca_rank, total_prop,fill = 0)  %>%
      gather(lca_rank, total_prop, -Query, -query_genus)  %>%
      group_by(Query) %>%
      mutate(lca_rank = factor(lca_rank, levels = c("species", "genus",
                                                    "family", "order",
                                                    "class","phylum",
                                                    "superkingdom"))) %>%
      arrange(lca_rank) %>%
      mutate(cum_prop = cumsum(total_prop))

#### Shigella EC Merge -------------------------------------------------------- 
shigella_class <- singleOrgMatchResults  %>%
      filter(query_genus == "Shigella", lca_name == "Enterobacteriaceae")  %>%
      .$match_taxid %>% unique()  %>% classification(db = 'ncbi')

e_coli_matches <- shigella_class  %>%
      map_df(bind_rows, .id = "query_id") %>%
      filter(name == "Escherichia coli")

singleOrgMatchResults_shigella_ecoli <- singleOrgMatchResults  %>%
      mutate(lca_name = if_else(match_taxid %in% e_coli_matches$query_id &
                                      query_genus == "Shigella",
                                "Echerichia coli", lca_name),
             lca_rank = if_else(match_taxid %in% e_coli_matches$query_id &
                                      query_genus == "Shigella",
                                "species", lca_rank),
             lca_id = if_else(match_taxid %in% e_coli_matches$query_id &
                                    query_genus == "Shigella",
                              "562", lca_id))

single_org_cum_ec <- singleOrgMatchResults_shigella_ecoli %>%
      filter(query_genus == "Shigella") %>%
      mutate(lca_rank = if_else(lca_rank %in% species_lvls,
                                "species", lca_rank)) %>%
      group_by(Query, lca_rank, query_genus) %>%
      summarise(total_prop = sum(`Final Guess`)) %>%
      spread(lca_rank, total_prop,fill = 0)  %>%
      gather(lca_rank, total_prop, -Query, -query_genus)  %>%
      group_by(Query) %>%
      mutate(lca_rank = factor(lca_rank, levels = c("species", "genus",
                                                    "family", "order",
                                                    "class","phylum",
                                                    "superkingdom"))) %>%
      arrange(lca_rank) %>%
      mutate(cum_prop = cumsum(total_prop))

## Species match lvl 
single_species_match <- single_org_cum %>% filter(lca_rank %in% c("species")) %>% 
      ungroup() %>% 
      mutate(query_genus = as.factor(query_genus), 
             query_genus_reorder = fct_reorder(query_genus, cum_prop, .desc = TRUE),
             match_thresh = if_else(cum_prop < 0.90,"darkorange","darkblue"))


## Species prop
spec_prop <- single_org_cum %>% filter(lca_rank %in% c("species")) %>% .$cum_prop
spec_0.999 <- sum(spec_prop <0.999)
spec_0.99 <- sum(spec_prop <0.99)
ec_st_shig_spec_prop <- single_org_cum %>% 
      filter(lca_rank %in% c("species"), 
             query_genus %in% c("Shigella","Escherichia","Staphlyococcus")) %>%
      .$cum_prop
ec_st_shig_spec_0.99 <- sum(ec_st_shig_spec_prop <0.99)
thresh_summary <- single_org_cum %>% filter(lca_rank == "species") %>% 
      mutate(thresh_0.99 = if_else(cum_prop < 0.99, 1,0),
             thresh_0.999 = if_else(cum_prop < 0.999, 1, 0)) %>% 
      group_by(query_genus) %>% 
      summarise(N = n(), below_0.99 = sum(thresh_0.99), below_0.999 = sum(thresh_0.999), 
                med_cum_prop = median(cum_prop))

non_ec_shig_st_med_spec <- single_org_cum %>% 
      filter(lca_rank %in% c("species"), 
             !query_genus %in% c("Shigella","Escherichia","Staphlyococcus")) %>%
      .$cum_prop %>% median()

shig_only <- single_org_cum %>% filter(query_genus == "Shigella") %>% 
      mutate(match_type = "shigella_only")

shig_and_esch <- single_org_cum_ec %>% mutate(match_type = "esch_and_shig")
shig_and_esch <- bind_rows(shig_only, shig_and_esch) %>% 
      mutate(query_group = paste(Query,match_type))
shig_match_change <- shig_and_esch %>% filter(lca_rank == "species") %>% 
      group_by(match_type) %>% 
      summarise(med_species = median(cum_prop)) %>% .$med_species

## Phage 
genus_add <- singleOrgMatchResults %>% select(Query, query_genus) %>% unique()
phage_total <- singleOrgResults %>% filter(grepl("phage",Genome)) %>% group_by(Query) %>% 
      summarise(total_phage_guess = sum(`Final Guess`)) %>% left_join(genus_add) %>% 
      ungroup() %>% 
      mutate(query_genus = fct_reorder(query_genus, total_phage_guess)) %>% 
      filter(query_genus != "Peptoclostridium")

# Individual contams -----------------------------------------------------------
clos_auto_top <- singleOrgResults %>% 
      filter(Query == "Clostridium_autoethanogenum_DSM_10061_uid219420") %>% top_n(1)
clos_auto_top_value <- clos_auto_top$`Final Guess`

ecoli_UMNK88 <- singleOrgResults %>% 
      filter(Query == "Escherichia_coli_UMNK88_uid42137")
ecoli_UMNK88_prov_value <- ecoli_UMNK88 %>% 
      filter(Genome == "ti|588|org|Providencia_stuartii") %>% .$`Final Guess`
ecoli_UMNK88_sal_value <- ecoli_UMNK88 %>% 
      filter(Genome == "ti|454169|org|Salmonella_enterica_subsp._enterica_serovar_Heidelberg_str._SL476" ) %>% 
            .$`Final Guess`

### Contam Table ---------------------------------------------------------------
contam_tbl_df <- contamSingleOrgMatchResults %>%
      mutate(lca_rank = if_else(lca_rank %in% species_lvls,
                                "species",lca_rank)) %>%
      filter(lca_rank %in% c("species","genus")) %>%
      group_by(Query, lca_rank) %>%
      summarize(total_prop = sum(`Final Guess`)) %>%
      spread(lca_rank,total_prop) %>%
      left_join(queryMeta) %>%
      mutate(`DNA length` = as.numeric(`DNA  length`),
             `DNA type` = if_else(`DNA length` < 1000000, "P","C")) %>%
      select(-GI, -Taxid) %>%
      group_by(Taxname, `DNA type`, species) %>%
      summarise(Acc = str_c(Accession, collapse = ", "),
                Mb = round(sum(`DNA length`)/1e6,digits = 2)) %>%
      gather("key","value", -Taxname, -`DNA type`, -species) %>%
      mutate(key = paste(`DNA type`, key)) %>%
      ungroup() %>% select(-`DNA type`) %>% spread(key,value) %>%
      rename(`Representative Strain` = Taxname, Species = species) %>%
      select(`Representative Strain`, Species ,
             `C Mb`, `C Acc`, `P Mb`, `P Acc`)


contamMixMatchResultsMin <- contamMixMatchResults %>%
      mutate(mix = as.numeric(mix), mix_contam = 1 - mix) %>%
      filter(lca_contam_rank %in%
                   c("genus","species","species group","subspecies")) %>%
      group_by(contam_name, target_name) %>% summarise(contam_min = min(mix_contam)) %>%
      mutate(contam_label = str_replace_all(contam_name, "_", " "),
             contam_label = str_replace(contam_label, " uid.*",""),
             contam_facet = str_sub(contam_name,1,4),
             target_facet = str_sub(target_name,1,4))

contam_min_df <- contamMixMatchResultsMin %>% ungroup() %>% 
      select(contam_label, target_facet, contam_min) %>% 
      spread(target_facet, contam_min)


contam_prop_df <- contamMixMatchResults %>% mutate(mix = as.numeric(mix), mix_contam = 1 - mix) %>%
      filter(lca_contam_rank %in%
                   c("genus","species","species group","subspecies")) %>%
      mutate(lca_contam_rank = ifelse(lca_contam_rank == "genus", "genus","species")) %>%
      group_by(contam_ds, target_name, contam_name, mix_contam) %>%
      summarise(contam_prop = sum(`Final Guess`)) %>%
      mutate(contam_facet = str_sub(contam_name,1,4),
             target_facet = str_sub(target_name,1,4),
             contam_label = str_replace_all(contam_name, "_", " "),
             contam_label = str_replace(contam_label, " uid.*","")) %>% 
      ungroup() %>% 
      select(contam_facet, target_facet, mix_contam, contam_prop)

contam_correlation <- contam_prop_df %>% 
      group_by(contam_facet, target_facet) %>% 
      summarise(prop_cor = cor(mix_contam, contam_prop))

contam_corr_quant <- contam_correlation$prop_cor %>% quantile(c(0.025,0.5, 0.975),na.rm = TRUE)

contam_resid <- contam_prop_df %>% 
      mutate(prop_resid = (contam_prop - mix_contam)/mix_contam) %>% 
      filter(mix_contam > 10^-5)

contam_resid_summary <- contam_resid %>% 
      group_by(contam_facet, target_facet) %>% 
      summarise(prop_resid_med = median(prop_resid), 
              prop_resid_min = min(prop_resid), 
              prop_resid_max = max(prop_resid),
              prop_resid_total = sum(prop_resid))

contam_resid_quant <- contam_resid_summary$prop_resid_total %>% 
      quantile(c(0.025,0.5, 0.975),na.rm = TRUE)

