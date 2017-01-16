## Genomic Purity Publication Variables, tables, and plots
library(xtable)
library(ProjectTemplate)

# load project
wd <- getwd()
setwd("../")
load.project()
setwd(wd)

##
## single genome data --------------------------------------------------------
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
      rename(Genus = query_genus, N = count) %>% arrange(Genus)


species_lvls <- c("species group", "species subgroup", "subspecies")
single_org_cum <- singleOrgMatchResults %>%
      filter(query_genus %in% c("Bacillus","Clostridium","Escherichia",
                                "Francisella","Listeria","Pseudomonas",
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
      .$match_taxid %>% unique() %>%
      {taxidClassification[names(taxidClassification) %in% .]}

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

## Check C. ljungdahlii also low
# singleOrgResults %>% filter(grepl("ljungdahlii", Query))


ecoli_UMNK88 <- singleOrgResults %>% 
      filter(Query == "Escherichia_coli_UMNK88_uid42137")
ecoli_UMNK88_prov_value <- ecoli_UMNK88 %>% 
      filter(Genome == "ti|588|org|Providencia_stuartii") %>% .$`Final Guess`
ecoli_UMNK88_sal_value <- ecoli_UMNK88 %>% 
      filter(Genome == "ti|454169|org|Salmonella_enterica_subsp._enterica_serovar_Heidelberg_str._SL476" ) %>% 
            .$`Final Guess`

### Contam Table ---------------------------------------------------------------
rep_strains <- c(
      "\\textit{Bacillus anthracis} str. Ames",                                       
      "\\textit{Clostridium botulinum} A str. Hall",                                  
      "\\textit{Escherichia coli} O157:H7 str. EC4115",                               
      "\\textit{Francisella tularensis} subsp. \\textit{tularensis} SCHU S4",                   
      "\\textit{Pseudomonas aeruginosa} PAO1",                                        
      "\\textit{Salmonella enterica} subsp. \\textit{enterica} serovar Typhimurium str. D23580",
      "\\textit{Staphylococcus aureus} subsp. \\textit{aureus} ED133",                          
      "\\textit{Yersinia pestis} CO92"
)

contam_tbl_raw_df <- contamSingleOrgMatchResults %>%
      mutate(lca_rank = if_else(lca_rank %in% species_lvls,
                                "species",lca_rank)) %>%
      filter(lca_rank %in% c("species","genus")) %>%
      group_by(Query, lca_rank) %>%
      summarize(total_prop = sum(`Final Guess`)) %>%
      spread(lca_rank,total_prop) %>%
      left_join(queryMeta) %>% 
      left_join(contamSingleOrgMapResults) %>% 
      mutate(`DNA length` = as.numeric(`DNA  length`),
             `DNA type` = if_else(`DNA length` < 1000000, "P","C"))

contam_tbl_df <- contam_tbl_raw_df%>%
      select(-GI, -Taxid, -`DNA type`) %>%
      group_by(Taxname, species, aligned_reads) %>%
      summarise(Mb = round(sum(`DNA length`)/1e6,digits = 2)) %>% 
      # gather("key","value", -Taxname, -species, -aligned_reads) %>%
      ungroup() %>% select(-Taxname) %>% add_column(Taxname = rep_strains) %>% 
      rename(`Representative Strain` = Taxname, 
             Species = species, 
             `Aligned Reads` = aligned_reads) %>%
      select(`Representative Strain`, Species , `Aligned Reads`,Mb)


contamMixMatchResultsMin <- contamMixMatchResults %>%
      mutate(mix = as.numeric(mix), mix_contam = 1 - mix) %>%
      filter(lca_contam_rank %in%
                   c("species","species group","subspecies")) %>%
      group_by(contam_name, target_name) %>% 
      summarise(contam_min = min(mix_contam)) %>%
      mutate(contam_label = str_replace_all(contam_name, "_", " "),
             contam_label = str_replace(contam_label, " uid.*",""),
             contam_facet = str_sub(contam_name,1,4),
             target_facet = str_sub(target_name,1,4))

contam_min_df <- contamMixMatchResultsMin %>% ungroup() %>% 
      select(contam_label, target_facet, contam_min) %>% 
      spread(target_facet, contam_min)


contam_prop_df <- contamMixMatchResults %>% 
      mutate(mix = as.numeric(mix), mix_contam = 1 - mix) %>%
      filter(lca_contam_rank %in%
                   c("species","species group","subspecies")) %>%
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

contam_corr_quant <- contam_correlation$prop_cor %>% 
      quantile(c(0.025,0.5, 0.975),na.rm = TRUE)

# ### Comparison of contaminant correlation in log and non-log space
# corr_comp <- contam_prop_df %>%
#       group_by(contam_facet, target_facet) %>%
#       summarise(logprop_cor = cor(log10(mix_contam + 1e-20), 
#                                   log10(contam_prop+ 1e-20))) %>%
#       left_join(contam_correlation) %>%
#       mutate(corr_diff = prop_cor - logprop_cor)
# 
# # Log space correlations slightly lower than non-log space correlations
# corr_comp$prop_cor %>% quantile(c(0.025,0.5, 0.975),na.rm = TRUE)
# corr_comp$logprop_cor %>% quantile(c(0.025,0.5, 0.975),na.rm = TRUE)
# corr_comp$corr_diff %>% quantile(c(0.025,0.5, 0.975),na.rm = TRUE)
# t.test(corr_comp$prop_cor, corr_comp$logprop_cor,
#        paired = TRUE,alternative = "greater")


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


## Contam Mix False Positives
contam_fp <- contamMixMatchResults %>% 
      filter(mix == "1.0", lca_contam_rank %in% c("species", species_lvls))
## No longer a false positive when looking for the contaminant at the species level
# contam_esch_df <- contam_fp %>% filter(target == "27739")
# contam_esch_prop <- contam_esch_df$`Final Guess`
# contam_esch_org <- "NEED TO CHECK" #taxidClassification[contam_esch_df$match_taxid][[1]]
# contam_esch_species <- "NEED TO CHECK" #contam_esch_org$name[9]

contam_yers_sal_df <- contam_fp %>% filter(target == "34", contaminant == "40625")
contam_yers_sal_prop <- contam_yers_sal_df$`Final Guess`
contam_yers_sal_org <- taxidClassification[contam_yers_sal_df$match_taxid][[1]]
contam_yers_sal_species <- contam_yers_sal_org$name[11]

contam_yers_esch_df <- contam_fp %>% filter(target == "34", contaminant == "27739")
contam_yers_esch_prop <- contam_yers_esch_df$`Final Guess`
contam_yers_esch_org <- taxidClassification[contam_yers_esch_df$match_taxid][[1]]
contam_yers_esch_species <- contam_yers_esch_org$name[length(contam_yers_esch_org$name)]


## Contam Single Org Read Counts
contam_single_read_count <- contamSingleOrgResults %>% group_by(Query) %>% 
      summarise(total_reads = sum(`Initial Best Hit Read Numbers`))


## Supplemental Table Baseline- List of genomes with metadata and accession numbers.
## Table description:: 
# Taxname - full organism name
# Taxid - NCBI taxonomic identifier
# GI - NCBI GenBank ID
# Accession - NCBI GenBank Accession number
# Length - size of DNA sequence in base pairs 
# Random - random number used to generate simulated sequence data

queryIdTbl %>% 
      left_join(queryMeta) %>% 
      left_join(singleOrgArt %>% rename(Query = org)) %>% 
      select(Taxname, Taxid, GI, Accession, `DNA  length`, rand) %>% 
      rename(Random = rand, Length = `DNA  length`) %>% 
      write_csv("supplemental_table_baseline_genomes.csv")

## Supplemental Table Contam - List of genomes with metadata and accession numbers
## Table description:: 
# Taxname - full organism name
# Taxid - NCBI taxonomic identifier
# GI - NCBI GenBank ID
# DNA type - C for chromosome and P for plasmid
# Accession - NCBI GenBank Accession number
# Length - size of DNA sequence in base pairs 
contam_tbl_raw_df %>% ungroup() %>% 
      select(Taxname, Taxid, GI, `DNA type`, Accession, `DNA length`) %>% 
      rename(Length = `DNA length`) %>% 
      write_csv("supplemental_table_contam_genomes.csv")
