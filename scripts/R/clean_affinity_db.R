###########################################################################
## This script processes tab-delimited Protein-Protein Interaction
## Affinity Database files into a more tidy format for easier processing 
## and annotation. Raw .tsv files are available here:
## 
## http://bmm.crick.ac.uk/~bmmadmin/Affinity
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## January 2016
###########################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

prots1 <- read_delim("affinity-1.tsv", "\t", skip = 1) # Affinity DB 1.0
prots2 <- read_delim("affinity-2.tsv", "\t", skip = 1) # Affinity DB 2.0

# Correct duplicate column names
names(prots1)[3] <- "pdb1"
names(prots1)[5] <- "pdb2"
names(prots2)[3] <- "pdb1"
names(prots2)[5] <- "pdb2"

# Select columns of interest
prots1 <- select(prots1, `Complex PDB`, Type, pdb1, pdb2)
prots2 <- select(prots2, `Complex PDB`, Type, pdb1, pdb2)

# Bind data frames togther and reshape
all_prots <- bind_rows(prots1, prots2) %>%
  gather(partner, partner_pdb, pdb1, pdb2) %>%
  select(-partner)

# Remove known enzyme complexes
all_prots <- filter(all_prots, str_detect(Type, "^O"))

# Export as CSV for manual annotation
write_csv(all_prots, "affinity_db_clean.csv")
