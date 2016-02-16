rm(list = ls())
library(cowplot)
library(dplyr)
library(readr)

calc_ROC <- function(data, radius) {
  false_neg <- filter(data, ACTIVE_SITE == 1, distance > radius) %>% nrow()
  true_pos <- filter(data, ACTIVE_SITE == 1, distance <= radius) %>% nrow()
  true_pos_rate <- true_pos/(true_pos + false_neg)
  
  true_neg <- filter(data, ACTIVE_SITE == 0, distance > radius) %>% nrow()
  false_pos <- filter(data, ACTIVE_SITE == 0, distance <= radius) %>% nrow()
  false_pos_rate <- false_pos/(false_pos + true_neg)
  
  data.frame(true_pos_rate, false_pos_rate)
}

master <- read_csv("master_data_table.csv", col_names=TRUE)

master.clean <- group_by(master[complete.cases(master),],pdb,Chain)
master.clean.wcn <- mutate(master.clean, distance = dist_wcn)

radii <- seq(0, 60, length.out=100)

roc1 <- lapply(radii, function(x) calc_ROC(master.clean, x)) %>% bind_rows() %>% mutate(model="distance")
roc2 <- lapply(radii, function(x) calc_ROC(master.clean.wcn, x)) %>% bind_rows() %>% mutate(model="wcn")
roc <- bind_rows(roc1, roc2)
roc_curve <- ggplot(roc, aes(x=false_pos_rate, y=true_pos_rate, color=factor(model, level=c("distance","wcn"), labels=c("Putative active site", "Maximum WCN")))) + 
  geom_line() +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  scale_color_discrete(name="Reference Point") +
  coord_fixed()
save_plot("roc_curve.pdf", roc_curve)
