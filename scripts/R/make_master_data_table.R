rm(list=ls())

library(readr)
library(dplyr)
library(tidyr)
library(mgcv)

calc_pred <- function(formula, data) {
  # Make model
  m <- lm(formula, data=data)
  # Pass predictions and residuals into dataframe
  ypred <- predict(m, data)
  # Bind new columns to data frame
  return(ypred)
}

calculate_dist <- function(file_name) {
  # Import protein CSV file and remove duplicate columns
  # Force column 'pdb' to character type
  protein_data <- read_csv(file_name, col_types = cols(pdb = "c"))
  protein_data <- protein_data[, unique(colnames(protein_data))] %>%
    rename(z_rate = rate, rate = raw_rate)

  # Tidy up our distance data by converting to a "long" table
  tidy_dist <- select(protein_data, Residue, Amino_Acid, ACTIVE_SITE, rate, wcnSC, RSA, matches("^[0-9]+[A-Z]$")) %>%
    gather(ref.site, distance, -rate, -rate.delta, -ACTIVE_SITE, -Residue, -Amino_Acid, -wcnSC, -RSA) %>%
    group_by(ref.site) 
  
  active_residues <- mutate(tidy_dist, ACTIVE_SITE=any(ACTIVE_SITE==1 & distance==0)) %>%
    filter(ACTIVE_SITE==TRUE)  
  
  if(nrow(active_residues) == 0) {
    active_residues <- data.frame(min_dist=rep(NA, nrow(protein_data)))
  } else {
    active_residues <- group_by(active_residues, Residue) %>% 
      summarize(min_dist=min(distance))
  }
  
  # Perform distance optimization to find best reference point
  tidy_dist <- mutate(tidy_dist, modelpred = calc_pred("rate ~ distance + wcnSC + RSA", data.frame(rate=rate, distance=distance, wcnSC=wcnSC, RSA=RSA)))
  
  max_cor <- summarize(tidy_dist, cor.dist = cor(distance, rate)) %>%
    filter(cor.dist == max(cor.dist))
  max_cor.k123 <- summarize(tidy_dist, cor.k123 = cor(distance, modelpred)) %>%
    filter(cor.k123 == max(cor.k123))

  # Pull out set of distances with highest correlation
  tidy_dist2 <- filter(tidy_dist, ref.site==max_cor$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_cor, by=c('ref.site'='ref.site'))
  tidy_dist.k123 <- filter(tidy_dist, ref.site==max_cor.k123$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_cor.k123, by=c('ref.site'='ref.site')) %>%
    rename(ref.site.k123 = ref.site, distance.k123 = distance)

  # Construct our final dataframe to output
  final_distances <- select(protein_data, -matches("^[0-9]+[A-Z]$")) %>%
      cbind(dist_active=active_residues$min_dist) %>%
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
             dist_active) %>%
        rename(pdb_site=SITE_NUMBER) %>%
        inner_join(tidy_dist2, by = c('Residue'='Residue')) %>%
        inner_join(tidy_dist.k123, by = c('Residue'='Residue'))
  
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
# test <- calculate_dist("merged/132L_merged.csv")
all_proteins <- lapply(paths, calculate_dist) %>% bind_rows()
write_csv(all_proteins, 'master_data_table.csv')
