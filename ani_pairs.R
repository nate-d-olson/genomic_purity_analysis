### Run dnadiff on sequence pairs

library(tidyverse)
library(stringr)

sim_orgs_dir <- list.dirs("data/simdir", recursive = FALSE)
db_dir <- "/media/nolson/genomic_purity/databases/Bacteria/"
sim_db_dir <- sim_orgs_dir %>% str_replace("data/simdir", db_dir) 
sim_orgs <- sim_orgs_dir %>% str_replace("data/simdir/","")
sim_dat <- list(org_name = sim_orgs, sim_dir = sim_orgs_dir, db_dir = sim_db_dir) %>% 
  transpose() %>% set_names(sim_orgs)

concat_fasta <- function(org_name, sim_dir, db_dir){
  fns <- list.files(db_dir, pattern = "*fna",full.names = TRUE)
  out_fasta <- paste0(sim_dir, "/", org_name, ".fna")
  if(!file.exists(out_fasta)) system2("cat",fns,stdout = out_fasta)
  out_fasta
}

fasta_files <- sim_dat %>% map(~concat_fasta(org_name = .$org_name, 
                                             sim_dir= .$sim_dir, 
                                             db_dir = .$db_dir))

# dir.create("data/pair_diff")
# 
# for(i in 1:(length(sim_orgs)-1)){
#   for(j in (i+1):length(sim_orgs)){
#     org1 <- sim_orgs[i]
#     org1_fna <- fasta_files[i]
#     org2 <- sim_orgs[j]
#     org2_fna <- fasta_files[j]
#     out_prefix <- paste0("data/pair_diff/", org1, "-", org2) 
#     system2("dnadiff", args = c("-p ", out_prefix, org1_fna, org2_fna))
#   }
# }

run_dnadiff <- function(idxs, sim_orgs, fasta_files){
  i <- idxs[1]
  org1 <- sim_orgs[i]
  org1_fna <- fasta_files[i]
  
  j <- idxs[2]
  org2 <- sim_orgs[j]
  org2_fna <- fasta_files[j]
  
  out_prefix <- paste0("data/pair_diff/", org1, "-", org2) 
  
  if(!file.exists(paste0(out_prefix,".report"))){
    system2("dnadiff", args = c("-p ", out_prefix, org1_fna, org2_fna))
    return("Run")
  }
  return("Not Run")
  
}

## Parallel implementation
library(parallel)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
genome_pairs <- combn(length(sim_orgs), 2, simplify = FALSE)
par_run_out <- parSapply(cl, genome_pairs, run_dnadiff, sim_orgs, fasta_files)
table(par_run_out)
