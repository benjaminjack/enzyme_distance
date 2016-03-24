rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

make_models_enzymes <- function(data) {
  lm1 <- lm(rate~wcnSC+RSA, data=data) %>% glance()
  lm2 <- lm(rate~wcnSC+RSA+dist_int, data=data) %>% glance()
  lm3 <- lm(rate~wcnSC+RSA+dist_active, data=data) %>% glance()
  data.frame(structure = lm1, struct_int = lm2, k123 = lm3)
}

make_models <- function(data) {
  lm1 <- lm(rate~wcnSC+RSA, data=data) %>% glance()
  lm2 <- lm(rate~wcnSC+RSA+dist_int, data=data) %>% glance()
  data.frame(structure = lm1, struct_int = lm2)
}

proteins <- read_csv("../master_data_table_pp.csv.gz")

proteins <- group_by(proteins, pdb, Chain)

shell_breaks_int <- c(0, seq(2.5, max(proteins$dist_int), by=5), max(proteins$dist_int))

# Append predictions from linear models and their respective residuals
proteins <- ungroup(proteins) %>%
  # Cut residues into shells as defined by shell_breaks above
  mutate(shell_int = cut(dist_int, labels=F, breaks = shell_breaks_int, include.lowest=T)) %>%
  mutate(shell_int = shell_int - 1) %>%
  group_by(pdb, Chain)
  
protein_mods <- proteins %>% do(mods = make_models(.)) %>% unnest() %>%
  select(pdb, Chain, structure.r.squared, struct_int.r.squared) %>%
  mutate(delta.int = struct_int.r.squared - structure.r.squared) %>%
  select(pdb, Chain, delta.int) %>%
  gather(comp, delta.r2, delta.int, -pdb, -Chain) %>%
  mutate(dataset="Non-enzymes")

enzymes <- read_csv("../master_data_table.csv.gz") %>%
  filter(multimer == 1) %>%
  mutate(interface = ifelse(RSA_mono - RSA > 0.10, T, F))

# Remove proteins without active site information
enzymes <- filter(enzymes, !is.na(dist_active), !is.na(dist_int)) %>% group_by(pdb,Chain)

shell_breaks <- c(0, seq(2.5, max(enzymes$dist_active), by=5), max(enzymes$dist_active))
shell_breaks_int <- c(0, seq(2.5, max(enzymes$dist_int), by=5), max(enzymes$dist_int))

# Append predictions from linear models and their respective residuals
enzymes <- ungroup(enzymes) %>%
  # Cut residues into shells as defined by shell_breaks above
  mutate(shell = cut(dist_active, labels=F, breaks = shell_breaks, include.lowest=T)) %>%
  mutate(shell = shell - 1) %>%
  mutate(shell_int = cut(dist_int, labels=F, breaks = shell_breaks_int, include.lowest=T)) %>%
  mutate(shell_int = shell_int - 1) %>%
  group_by(pdb, Chain)

enzymes_mods <- enzymes %>% do(mods = make_models_enzymes(.)) %>% unnest() %>%
  select(pdb, Chain, structure.r.squared, struct_int.r.squared, k123.r.squared) %>%
  mutate(delta.int = struct_int.r.squared - structure.r.squared, delta.act = k123.r.squared - structure.r.squared) %>%
  select(pdb, Chain, delta.int, delta.act) %>%
  rename(delta.int.enzyme = delta.int) %>%
  gather(comp, delta.r2, delta.int.enzyme:delta.act, -pdb, -Chain) %>%
  mutate(dataset="Enzymes")

cats <- filter(enzymes, ACTIVE_SITE == T) %>%
  mutate(type="cat", dataset="enzyme") %>%
  select(type, rate, dataset)
ints <- filter(enzymes, interface == T) %>%
  mutate(type="interface_enzyme", dataset="enzyme") %>%
  select(type, rate, dataset)
enzyme_nonint <- filter(enzymes, interface == F, ACTIVE_SITE == F) %>%
  mutate(type="no_cat_int", dataset="enzyme") %>%
  select(type, rate, dataset)
proteins_int <- filter(proteins, interface == T) %>%
  mutate(type="interface_protein", dataset="prot") %>%
  select(type, rate, dataset)
proteins_noint <- filter(proteins, interface == F) %>%
  mutate(type="no_interface_protein", dataset="prot") %>%
  select(type, rate, dataset)

cats_ints <- bind_rows(cats, ints, enzyme_nonint, proteins_int, proteins_noint) %>%
  mutate(dataset = factor(dataset, levels=c("enzyme", "prot"), labels=c("Enzymes", "Non-enzymes")))

plot1 <- ggplot(cats_ints, aes(x=factor(type, 
                               levels=c("cat", "interface_enzyme", "no_cat_int", "interface_protein", "no_interface_protein"),
                               labels=c("Catalytic",
                                        "Interface",
                                        "All other sites",
                                        "Interface",
                                        "All other sites")), 
                               fill=dataset, y=rate)) +
  coord_cartesian(ylim=c(0,1.1)) +
  stat_summary(fun.y=mean, geom="bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  xlab("Residue Type") +
  ylab("Relative Rate") +
  scale_fill_grey(name="Data set") + 
  facet_grid(.~dataset, scale="free_x", space="free_x") +
  theme(legend.position="none")

# Write data to file for plot1
data_out <- cats_ints %>% group_by(type, dataset) %>% summarize(mean=mean(rate), se=sd(rate)/sqrt(n()))
write_csv(data_out, "fig_enzyme_a.csv")

enzyme_rate <- ggplot(filter(enzymes, shell_int <= 10), aes(x=factor(shell_int), y=rate)) + 
   geom_violin(scale="width", trim=TRUE, aes(fill=..count..)) +
   stat_summary(fun.y=mean, geom="point") +
   stat_summary(fun.y=mean, geom="line", aes(group=1)) +
   coord_cartesian(ylim=c(0,3)) +
   xlab('Shell') +
   ylab('Relative Rate') +
   scale_fill_gradient(name="Residue\nCount", low="lightblue", high="blue", trans="sqrt") +
   theme(legend.key.height = unit(25, "pt"))

non_enzyme_rate <- ggplot(filter(proteins, shell_int <= 10), aes(x=factor(shell_int), y=rate)) + 
  geom_violin(scale="width", trim=TRUE, aes(fill=..count..)) +
  stat_summary(fun.y=mean, geom="point") +
  stat_summary(fun.y=mean, geom="line", aes(group=1)) +
  coord_cartesian(ylim=c(0,3)) +
  xlab('Shell') +
  ylab('Relative Rate') +
  scale_fill_gradient(name="Residue\nCount", low="lightblue", high="blue", trans="sqrt") +
  theme(legend.key.height = unit(25, "pt"))

# save_plot("./fig_nonenzyme.pdf", plot1, base_width = 8, base_height = 5)

all_data <- bind_rows(enzymes_mods, protein_mods)

data_out <- all_data %>% group_by(comp, dataset) %>% summarize(mean=mean(delta.r2), se=sd(delta.r2)/sqrt(n()))
write_csv(data_out, "fig_enzyme_b_data.csv")

delta_r2_plot <- ggplot(all_data, aes(x=factor(comp, 
                              levels = c('delta.act', 
                                         'delta.int.enzyme', 
                                         'delta.int'),
                              labels = c('Catalytic',
                                         'Interface\n(enzymes)',
                                         'Interface\n(non-enzymes)')), y=delta.r2)) + 
  stat_summary(fun.y=mean, geom="bar", aes(fill=dataset)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_path(aes(x=x, y=y), data=data.frame(x=c(1,1,3,3), y=c(0.060, 0.062, 0.062, 0.060))) +
  geom_text(aes(x=2, y=0.0624, label="***"), size=7) +
  geom_path(aes(x=x, y=y), data=data.frame(x=c(1,1,2,2), y=c(0.056, 0.058, 0.058, 0.056))) +
  geom_text(aes(x=1.5, y=0.0584, label="***"), size=7) +
  ylab(expression(paste(Delta, R^{"2"}))) +
  scale_fill_grey() +
  theme(legend.position="none") +
  xlab("Distance Model")

fig1 <- ggdraw() +
  draw_plot(plot1, x=0, y=0.5, width=0.6, height=0.5) +
  draw_plot(delta_r2_plot, x=0.6, y=0.475, width=0.4, height=0.5) +
  draw_plot(enzyme_rate, x=0, y=0, width=0.5, height=0.47) +
  draw_plot(non_enzyme_rate, x=0.5, y=0, width=0.5, height=0.47) +
  draw_plot_label(c('a','b','c','d'), x = c(0, 0.6, 0, 0.5), y= c(1,1,0.47,0.47))

# plot_grid(plot1, delta_r2_plot, enzyme_rate, non_enzyme_rate, align="h")

save_plot("fig_nonenzyme.pdf", fig1, base_width = 10, base_height = 7.5)

all_data %>% group_by(comp, dataset) %>% summarize(mean=mean(delta.r2))

t.test(filter(all_data, comp == "delta.act")$delta.r2, filter(all_data, comp == "delta.int.enzyme")$delta.r2)
t.test(filter(all_data, comp == "delta.act")$delta.r2, filter(all_data, comp == "delta.int")$delta.r2)

cats_ints %>% group_by(type, dataset) %>% summarize(mean=mean(rate))

t.test(filter(cats_ints, type == "cat")$rate, filter(cats_ints, type == "interface_enzyme")$rate)
t.test(filter(cats_ints, type == "cat")$rate, filter(cats_ints, type == "no_cat_int")$rate)
t.test(filter(cats_ints, type == "interface_enzyme")$rate, filter(cats_ints, type == "no_cat_int")$rate)
t.test(filter(cats_ints, type == "interface_protein")$rate, filter(cats_ints, type == "no_interface_protein")$rate)
