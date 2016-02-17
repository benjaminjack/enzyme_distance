library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(tidyr)

enzymes <- read_csv("./data/enzymes/master_data_table.csv")
prots <- read_csv("./data/nonenzymes/master_data_table_pp.csv")

enzymes %>% group_by(pdb, Chain) %>%
  do(mod.structure = glance(lm(rate ~ RSA + wcnSC, data = .)), r.k123 = unique(.$cor.k123)) %>%
  mutate(r.k123 = r.k123[[1]]^2) %>%
  unnest() %>%
  select(pdb, Chain, r.squared, r.k123) %>% 
  ungroup() %>%
  summarize(mean_r = mean(r.squared), mean_r.k123 = mean(r.k123)) %>%
  gather(var, mean) %>%
  mutate(dataset = "enzymes") -> enzymes2

prots %>% group_by(pdb, Chain) %>%
  do(mod.structure = glance(lm(rate ~ RSA + wcnSC, data = .)), r.k123 = unique(.$cor.k123)) %>%
  mutate(r.k123 = r.k123[[1]]^2) %>%
  unnest() %>%
  select(pdb, Chain, r.squared, r.k123) %>%
  ungroup() %>%
  summarize(mean_r = mean(r.squared), mean_r.k123 = mean(r.k123)) %>%
  gather(var, mean) %>%
  mutate(dataset = "proteins") -> prots2

all_data <- bind_rows(enzymes2, prots2)

ggplot(all_data, aes(x=dataset, y=mean, fill=factor(var, levels=c("mean_r", "mean_r.k123"), labels=c("structure", "structure + d")))) + 
  geom_bar(stat="identity", position="dodge") +
  ylab("mean r^2") +
  scale_fill_discrete(name="model")
  