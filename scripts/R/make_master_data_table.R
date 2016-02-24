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

calculate_dist <- function(file_name) {
  # Import protein CSV file and remove duplicate columns
  # Force column 'pdb' to character type
  protein_data <- read_csv(file_name, col_types = cols(pdb = "c"))
  protein_data <- protein_data[, unique(colnames(protein_data))] %>%
    rename(z_rate = rate) %>% 
    mutate(max_wcn = ifelse(wcnSC == max(wcnSC), 1, 0)) %>%
    mutate(rate = raw_rate/mean(raw_rate))
  protein_data <- mutate(protein_data, interface = ifelse(RSA_mono - RSA > 0.10, T, F))
  
  # Tidy up our distance data by converting to a "long" table
  tidy_dist <- select(protein_data, Residue, Amino_Acid, interface, ACTIVE_SITE, max_wcn, rate, wcnSC, RSA, matches("^[0-9]+[A-Z]$")) %>%
    gather(ref.site, distance, -rate, -ACTIVE_SITE, -max_wcn, -interface, -Residue, -Amino_Acid, -wcnSC, -RSA) %>%
    group_by(ref.site) 
  
  active_residues <- mutate(tidy_dist, ACTIVE_SITE=any(ACTIVE_SITE==1 & distance==0)) %>%
    filter(ACTIVE_SITE==TRUE)
  
  interface_res <- mutate(tidy_dist, interface=any(interface==1 & distance==0)) %>%
    filter(interface==TRUE)
  
  max_wcn_residues <- mutate(tidy_dist, max_wcn=any(max_wcn==1 & distance==0)) %>%
    filter(max_wcn==TRUE) %>% 
    rename(ref.site.wcn = ref.site, dist_wcn = distance) %>%
    select(Residue, ref.site.wcn, dist_wcn)
  
  if(nrow(active_residues) == 0) {
    active_residues <- data.frame(min_dist=rep(NA, nrow(protein_data)))
  } else {
    active_residues <- group_by(active_residues, Residue) %>% 
      summarize(min_dist=min(distance))
  }
  
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

  # Pull out set of distances with highest correlation
  tidy_dist2 <- filter(tidy_dist, ref.site==max_cor$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_cor, by=c('ref.site'='ref.site'))
  tidy_dist.k123 <- filter(tidy_dist, ref.site==max_r2_k123$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_r2_k123, by=c('ref.site'='ref.site')) %>%
    rename(ref.site.k123 = ref.site, distance.k123 = distance)

  # Construct our final dataframe to output
  final_distances <- select(protein_data, -matches("^[0-9]+[A-Z]$")) %>%
      cbind(dist_active=active_residues$min_dist) %>%
      cbind(dist_int=interface_res$min_dist) %>%
      select(pdb,
             Chain,
             multimer,
             Residue, 
             SITE_NUMBER, 
             Amino_Acid,
             rate,
             z_rate,
             RSA,
             RSA_mono,
             wcnSC,
             wcnCA,
             wcnSC_mono,
             wcnCA_mono,
             ACTIVE_SITE,
             dist_active,
             interface,
             dist_int) %>%
        rename(pdb_site=SITE_NUMBER) %>%
        inner_join(tidy_dist2, by = c('Residue'='Residue')) %>%
        inner_join(tidy_dist.k123, by = c('Residue'='Residue')) %>% 
        inner_join(max_wcn_residues, by = c('Residue'='Residue'))
  
  # Make model
  if (!is.na(final_distances$dist_active)) {
    m2 <- gam(rate ~ s(dist_active, bs = "cs"), data=final_distances)
    preds <- predict(m2, data.frame(dist_active=final_distances$dist_active))
  } else {
    preds <- rep(NA, nrow(final_distances))
  }
  final_distances <- cbind(final_distances, data.frame(gam_rate=preds))
  
  print(file_name)
  return(final_distances)
}

paths <- dir("merged", "\\.csv$", full.names = TRUE)
# test <- calculate_dist("merged/1FR2_merged.csv")
all_proteins <- lapply(paths, calculate_dist) %>% bind_rows()
write_csv(all_proteins, '../../master_data_table.csv')
