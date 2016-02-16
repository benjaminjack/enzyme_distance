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

load_structure_props <- function(file_name) {
  # Import protein CSV file and remove duplicate columns
  protein_data <- read_csv(file_name)
  protein_data <- protein_data[, unique(colnames(protein_data))] %>%
    rename(z_rate = rate) %>%
    mutate(rate = raw_rate/mean(raw_rate))
  
  pdb_name <- data.frame(pdb = rep(str_match(basename(file_name), "([A-Z0-9]{4})_")[2], nrow(protein_data)))
  
  protein_data <- bind_cols(protein_data, pdb_name)
  # return(protein_data)
  # Tidy up our distance data by converting to a "long" table
  tidy_dist <- select(protein_data, Residue, Amino_Acid, rate, wcnSC, RSA, matches("^[0-9]+[A-Z]$")) %>%
    gather(ref.site, distance, -rate, -Residue, -Amino_Acid, -wcnSC, -RSA) %>%
    group_by(ref.site)
  
  # Perform distance optimization to find best reference point
  tidy_dist <- mutate(tidy_dist, modelpred = calc_pred("rate ~ distance + RSA + wcnSC", data.frame(rate=rate, distance=distance, wcnSC=wcnSC, RSA=RSA)))
  
  # return(tidy_dist)
  
  max_cor.k123 <- summarize(tidy_dist, cor.k123 = cor(rate, modelpred)) %>%
    filter(cor.k123 == max(cor.k123))
  
  tidy_dist.k123 <- filter(tidy_dist, ref.site==max_cor.k123$ref.site) %>%
    select(Residue, ref.site, distance) %>%
    inner_join(max_cor.k123, by=c('ref.site'='ref.site')) %>%
    rename(ref.site.k123 = ref.site, distance.k123 = distance)
  
  # Construct our final dataframe to output
  final_distances <- select(protein_data, -matches("^[0-9]+[A-Z]$")) %>%
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
             wcnCA_mono) %>%
      inner_join(tidy_dist.k123, by = c('Residue'='Residue'))
  
  print(file_name)
  return(final_distances)
}

paths <- dir("./mapped", "\\.csv$", full.names = TRUE)
# test <- load_structure_props("mapped/1A2K_A_merged.csv")
all_proteins <- lapply(paths, load_structure_props) %>% bind_rows()
write_csv(all_proteins, 'master_data_table_pp.csv')
