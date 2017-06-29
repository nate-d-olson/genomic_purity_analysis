### Removing duplicate baseline genome entries 
supplemental_table_baseline_genomes <- read_csv("../reports/supplemental_table_baseline_genomes.csv")
duplicate_entries <- supplemental_table_baseline_genomes %>% 
      group_by(Taxname, Taxid, GI,Accession,Length) %>% 
      summarise(count = n()) %>% filter(count > 1)

## Removed extra log files using 
# Bacillus files using 
# rm Bacillus_*/logs/art_sim-2015-06-19-18*
# rm Bacillus_*/logs/pathomap-2015-06-19-18*
# E. coli 
# rm Escherichia_coli_*/logs/art_sim-2015-06-09*log
# rm Escherichia_coli_*/logs/art_sim-2015-06-10-11-18*log
## Missing a log file for Escherichia coli ....
### 62 pathoid log files now only 61 art_sim log files 
## Copied file back from copy of files from repo 
## Salmonella
# rm Salmonella_*/logs/art_sim-2015-07-04-00-24*log
# rm Salmonella_*/logs/art_sim-2015-07-04-00-25*log
# rm Salmonella_*/logs/art_sim-2015-07-04-00-26*log
# rm Salmonella_enterica_arizonae_serovar_62_z4_z23__RSK2980_uid13030/logs/art_sim-2015-07-04-00-32-36.log

