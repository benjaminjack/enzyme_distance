proteins <- read_csv("./data/nonenzymes/master_data_table_pp.csv") %>%
  mutate(interface = ifelse(RSA_mono - RSA > 0.10, T, F))

enzymes <- read_csv("./data/enzymes/master_data_table.csv") %>%
  mutate(interface = ifelse(RSA_mono - RSA > 0.10, T, F))

cats <- filter(enzymes, ACTIVE_SITE == T) %>%
  mutate(type="cat", dataset="enzyme") %>%
  select(type, rate, dataset)
ints <- filter(enzymes, interface == T) %>%
  mutate(type="interface_enzyme", dataset="enzyme") %>%
  select(type, rate, dataset)
enzyme_nonint <- filter(enzymes, interface == F, ACTIVE_SITE == F) %>%
  mutate(type="no_cat_int", dataset="enzyme") %>%
  select(type, rate, dataset)
proteins_int <- filter(enzymes, interface == T) %>%
  mutate(type="interface_protein", dataset="prot") %>%
  select(type, rate, dataset)
proteins_noint <- filter(enzymes, interface == F) %>%
  mutate(type="no_interface_protein", dataset="prot") %>%
  select(type, rate, dataset)


cats_ints <- bind_rows(cats, ints, enzyme_nonint, proteins_int, proteins_noint)

plot1 <- ggplot(cats_ints, aes(x=factor(type, 
                               levels=c("cat", "interface_enzyme", "no_cat_int", "interface_protein", "no_interface_protein"),
                               labels=c("Catalytic",
                                        "Interface",
                                        "All other sites",
                                        "Interface\n(non-enzyme)",
                                        "All other sites\n(non-enzyme)")), fill=dataset, y=rate)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_summary(fun.y=mean, geom="bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  xlab("Residue Type") +
  ylab("Relative Rate")

save_plot("./figures/fig_nonenzyme.pdf", plot1.2, ncol = 2)


