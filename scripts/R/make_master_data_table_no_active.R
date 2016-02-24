rm(list=ls())

library(readr)
library(dplyr)
library(tidyr)
library(mgcv)
library(ic.infer)
library(broom)

make_model <- function(data) {
  # Make model
  m <- lm(rate ~ distance + RSA + wcnSC, data=data)
  # Pass predictions and residuals into dataframe
  ui <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, -1))
  ci <- c(0, 0, 0)
  orm <- orlm(m, ui=ui, ci=ci)
  # Bind new columns to data frame
  data.frame(r2.k123=orm$R2, dist.coeff=orm$b.restr["distance"])
}

load_structure_props <- function(file_name) {
  # Import protein CSV file and remove duplicate columns
  protein_data <- read_csv(file_name)
  protein_data <- protein_data[, unique(colnames(protein_data))] %>%
    rename(z_rate = rate) %>%
    mutate(rate = raw_rate/mean(raw_rate))
  protein_data <- mutate(protein_data, interface = ifelse(RSA_mono - RSA > 0.10, T, F))
  pdb_name <- data.frame(pdb = rep(str_match(basename(file_name), "([A-Z0-9]{4})_")[2], nrow(protein_data)))
  
  protein_data <- bind_cols(protein_data, pdb_name)
  # return(protein_data)
  # Tidy up our distance data by converting to a "long" table
  tidy_dist <- select(protein_data, Residue, Amino_Acid, interface, rate, wcnSC, RSA, matches("^[0-9]+[A-Z]$")) %>%
    gather(ref.site, distance, -rate, -Residue, -Amino_Acid, -interface, -wcnSC, -RSA) %>%
    group_by(ref.site)
  
  interface_res <- mutate(tidy_dist, interface=any(interface==1 & distance==0)) %>%
    filter(interface==TRUE)
  
  if(nrow(interface_res) == 0) {
    interface_res <- data.frame(min_dist=rep(NA, nrow(protein_data)))
  } else {
    interface_res <- group_by(interface_res, Residue) %>% 
      summarize(min_dist=min(distance))
  }
  
  # Perform distance optimization to find best reference point
  max_r2_k123 <- do(tidy_dist, make_model(.)) %>% ungroup() %>% filter(r2.k123 == max(r2.k123))
  # return(tidy_dist)
  
  max_cor <- summarize(tidy_dist, r.dist = cor(distance, rate)^2) %>%
    filter(r.dist == max(r.dist))
  
  tidy_dist2 <- filter(tidy_dist, ref.site==max_cor$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_cor, by=c('ref.site'='ref.site'))
  
  tidy_dist.k123 <- filter(tidy_dist, ref.site==max_r2_k123$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_r2_k123, by=c('ref.site'='ref.site')) %>%
    rename(ref.site.k123 = ref.site, distance.k123 = distance)
  
  # Construct our final dataframe to output
  final_distances <- select(protein_data, -matches("^[0-9]+[A-Z]$")) %>%
      cbind(dist_int=interface_res$min_dist) %>%
      select(pdb,
             Chain,
             Residue, 
             Amino_Acid,
             rate,
             z_rate,
             RSA,
             RSA_mono,
             wcnSC,
             wcnCA,
             wcnSC_mono,
             wcnCA_mono,
             interface,
             dist_int) %>%
      inner_join(tidy_dist.k123, by = c('Residue'='Residue')) %>%
      inner_join(tidy_dist2, by = c('Residue'='Residue'))
  
  print(file_name)
  return(final_distances)
}

paths <- dir("./mapped", "\\.csv$", full.names = TRUE)
# test <- load_structure_props("mapped/1EER_A_merged.csv")
all_proteins <- lapply(paths, load_structure_props) %>% bind_rows()
write_csv(all_proteins, '../../master_data_table_pp.csv')
