rm(list = ls())

library(dplyr)
library(tidyr)
library(broom)

make_models <- function(data) {
  
  # Make all of our models, "tidy" them up, and label with the formula
  models <- list(
    tidy(lm(rate~RSA, data)) %>% mutate(formula = "rate~RSA"),
    tidy(lm(rate~wcnSC, data)) %>% mutate(formula = "rate~wcnSC"),
    tidy(lm(rate~wcnSC+RSA, data)) %>% mutate(formula = "rate~wcnSC+RSA"),
    tidy(lm(rate~dist_active, data)) %>% mutate(formula = "rate~dist_active"),
    tidy(lm(rate~wcnSC+RSA+dist_active, data)) %>% mutate(formula = "rate~wcnSC+RSA+dist_active"),
    tidy(lm(rate~wcnSC+dist_active, data)) %>% mutate(formula = "rate~wcnSC+dist_active"),
    tidy(lm(rate~RSA+dist_active, data)) %>% mutate(formula = "rate~RSA+dist_active"),
    tidy(lm(rate~RSA_mono, data)) %>% mutate(formula = "rate~RSA_mono"),
    tidy(lm(rate~wcnSC_mono, data)) %>% mutate(formula = "rate~wcnSC_mono"),
    tidy(lm(rate~wcnSC_mono+RSA_mono, data)) %>% mutate(formula = "rate~wcnSC_mono+RSA_mono"),
    tidy(lm(rate~wcnSC_mono+RSA_mono+dist_active, data)) %>% mutate(formula = "rate~wcnSC_mono+RSA_mono+dist_active"),
    tidy(lm(rate~wcnSC_mono+dist_active, data)) %>% mutate(formula = "rate~wcnSC_mono+dist_active"),
    tidy(lm(rate~RSA_mono+dist_active, data)) %>% mutate(formula = "rate~RSA_mono+dist_active")
  )
  
  bind_rows(models)
  
}

master <- read_csv("master_data_table.csv")

master %>% filter(!is.na(dist_active)) %>% # Only include proteins with active site info
  group_by(pdb, Chain) %>%
  do(mod = make_models(.)) %>%
  unnest() %>%
  select(pdb, Chain, formula, term, estimate, std.error, statistic, p.value) -> model_table

write_csv(model_table, "model_parameters.csv")


