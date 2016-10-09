## Genomic Purity Publication Variables, tables, and plots
library(knitr)
library(xtable)
library(forcats)
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
thresh_summary <- single_org_cum %>% filter(lca_rank == "species") %>% mutate(thresh_0.99 = if_else(cum_prop < 0.99, 1,0),
                                                                              thresh_0.999 = if_else(cum_prop < 0.999, 1, 0)) %>% 
      group_by(query_genus) %>% 
      summarise(N = n(), below_0.99 = sum(thresh_0.99), below_0.999 = sum(thresh_0.999), med_cum_prop = median(cum_prop))

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

## Contam Table ---------------------------------------------------------------
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
contam_min_df<- contamMixMatchResultsMin %>% ungroup() %>% 
      select(contam_label, target_facet, contam_min) %>% 
      spread(target_facet, contam_min)

contam_plot <-contamMixMatchResultsMin %>% ggplot(aes(x = mix_contam, y = contam_prop)) +
      geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "grey60") +
      geom_vline(data = contamMixMatchResultsMin,
                 aes(xintercept = contam_min, color = contam_label), linetype = 3) +
      geom_path(aes(color = contam_label), alpha = 0.5) +
      geom_point(aes(color = contam_label)) +
      facet_grid(contam_facet~target_facet) + scale_x_log10() + scale_y_log10() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
      guides(color = guide_legend(ncol = 2)) +
      labs(x = "Contaminant Proportion",
           y = "Proportion of Matched Reads",
           color = "Contaminant")