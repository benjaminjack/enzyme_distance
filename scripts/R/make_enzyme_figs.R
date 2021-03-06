rm(list=ls())
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(cowplot)

# =============================================================
# ANALYSIS OPTIONS
# =============================================================

# These are for testing purposes only. These options are now all passed
# as arguments in a function.

# DATA_TABLE <- "master_data_table.csv" # Path to master data table
# RAW_RATES <- T # Use raw rate4site rates rather than z-scores
# MULTIMERS_ONLY <- F # Use only multimeric proteins
# MONOMERS_ONLY <- F  # Use only monomeric proteins
# MULTI_RSA <- T # Use multimeric RSA
# MULTI_WCN <- T # Use multimeric WCN
# NO_GLY <- F         # Remove glycines before analysis
# NO_INTERFACE <- F   # Remove interface residues

# =============================================================
# END OPTIONS
# =============================================================

# Set some defaults for ggplot2
theme_set(
  theme_cowplot() + 
    theme(strip.background = element_rect(colour="black", 
                                          fill="grey80", 
                                          size=0.5, 
                                          linetype = "solid"))
)

# =============================================================
# DEFINE FUNCTIONS FOR ANALYSIS
# =============================================================

# Function to calculate correlations given a linear model or just two variables
calc_cor <- function(formula=NA, data) {
  if(is.na(formula)) {
    # No formula given, so calculate pearson correlation without linear model
    return(cor(data$y, data$x1))
  }
  m <- lm(formula, data=data)
  ypred <- predict(m, data)
  return(cor(data$y, ypred))
}

# Determine coefficient from a linear model
calc_coeff <- function(formula, data) {
  m <- lm(formula, data=data)
  return(coef(m)[2])
}

calc_pred <- function(formula, combo, data) {
  # Make model
  m <- lm(formula, data=data)
  # Pass predictions and residuals into dataframe
  ypred <- data.frame(predict(m, data), residuals(m))
  # Make column names
  title <- paste(combo, collapse=".")
  colnames(ypred) <- c(str_c(title, ".pred"), str_c(title, ".res"))
  # Bind new columns to data frame
  return(bind_cols(data,ypred))
}

make_lms <- function(df, y, vars) {
  # Initiate list to hold combinations of variables
  x <- c()
  # Create all combos of variables for linear models
  for (n in seq(length(vars))) { x <- append(x,combn(vars,n,simplify = F)) }
  # Loop over combos
  for (combo in x) {
    # Generate formula for linear model
    formula <- as.formula(paste(y,'~',paste(combo, collapse="+")))
    # Compute linear model predictions by protein
    df <- do(df, calc_pred(formula, combo, .))
  }
  return(df)
}

# Assign null model groups for facetting
get_null_group <- function(pred) {
  pred <- as.character(pred)
  pred[pred == 'wcnSC.pred'] <- 2
  pred[pred == 'wcnSC.dist_active.pred'] <- 2
  pred[pred == 'wcnSC.res'] <- 2
  pred[pred == 'wcnSC.dist_active.res'] <- 2
  pred[pred == 'RSA.pred'] <- 1
  pred[pred == 'RSA.dist_active.pred'] <- 1
  pred[pred == 'RSA.res'] <- 1
  pred[pred == 'RSA.dist_active.res'] <- 1
  pred[pred == 'rate'] <- 0
  pred[pred == 'rate1'] <- 1
  pred[pred == 'rate2'] <- 2
  pred[pred == 'rate3'] <- 3
  pred[pred == 'dist_active.pred'] <- 1
  pred[pred == 'dist_active.pred2'] <- 2
  pred[pred == 'dist_active.pred3'] <- 3
  pred[pred == 'dist_active.res'] <- 1
  pred[pred == 'dist_active.res2'] <- 2
  pred[pred == 'dist_active.res3'] <- 3
  pred[pred == 'wcnSC.RSA.pred'] <- 3
  pred[pred == 'wcnSC.RSA.dist_active.pred'] <- 3
  pred[pred == 'wcnSC.RSA.res'] <- 3
  pred[pred == 'wcnSC.RSA.dist_active.res'] <- 3
  return(pred)
}

# Assign models as null or alternate
get_null <- function(pred) {
  pred <- as.character(pred)
  pred <- str_replace(pred, '^dist_active.*', 'd')
  pred <- str_replace(pred, '.*rate.*', 'rate')
  pred <- str_replace(pred, '.*\\.dist_active.*', 'alt')
  pred[pred != 'alt' & pred != 'd' & pred != 'rate'] <- 'null'
  return(pred)
}

# For labelling facets
label_size_group <- as_labeller(c(`1` = "Small Proteins",
                                  `2` = "Medium Proteins",
                                  `3` = "Large Proteins"))
label_null_group <- as_labeller(c(`1` = "RSA",
                                  `2` = "WCN",
                                  `3` = "WCN + RSA"))
label_rsa_group <- as_labeller(c(`1` = "Core",
                                 `2` = "Intermediate",
                                 `3` = "Surface"))

# 
# ALL OF THE FIGURE MAKING CODE IS WRAPPED IN THE FOLLOWING FUNCTION:
#
 
make_figures <- function(master.clean, # Path to master data table
                         RAW_RATES = T, # Use raw rate4site rates rather than z-scores
                         MULTIMERS_ONLY = F, # Use only multimeric proteins
                         MONOMERS_ONLY = F,  # Use only monomeric proteins
                         MULTI_RSA = T, # Use multimeric RSA
                         MULTI_WCN = T, # Use multimeric WCN
                         NO_GLY = F,         # Remove glycines before analysis
                         NO_INTERFACE = F,   # Remove interface residues
                         suffix = "") {
  
  # =============================================================
  # DATA CLEANING
  # =============================================================
  
  # Use raw rates ?
  if (RAW_RATES == F) {
    master.clean <- mutate(master.clean, rate = z_rate)
  }
  
  # Interface residues ?
  if (NO_INTERFACE == T) {
    master.clean <- filter(master.clean, abs(RSA_mono - RSA) < 0.10)
  }
  
  # Multimers only ?
  if (MULTIMERS_ONLY == T) {
    master.clean <- filter(master.clean, multimer == 1)
  }
  
  # Monomers only ?
  if (MONOMERS_ONLY == T) {
    master.clean <- filter(master.clean, multimer == 0)
  }
  
  # Multimeric WCN ?
  if (MULTI_WCN == F) {
    master.clean <- mutate(master.clean, wcnSC = wcnSC_mono)
  }
  
  # Multimeric RSA ?
  if (MULTI_RSA == F) {
    master.clean <- mutate(master.clean, RSA = RSA_mono)
  }
  
  # Glycines ?
  if (NO_GLY == T) {
    master.clean <- filter(master.clean, Residue != 'G')
  }
  
  # =============================================================
  # DATA PROCESSING
  # =============================================================
  
  # Append predictions from linear models and their respective residuals
  master.clean <- make_lms(master.clean, "rate", c("wcnSC", "RSA", "dist_active"))
  
  master.lmcor <- group_by(master.clean, pdb, multimer) %>% 
    summarise_each(funs(cor = cor.test(., rate)$estimate), ends_with(".pred"))
  
  # =============================================================
  # FIGURE 1: RATE DISTRIBUTIONS
  # =============================================================
  
  # Basic distance vs rate plots
  p1a <- ggplot(master.clean, aes(x=dist_active, y=rate)) + 
    ylab('Relative Rate') + 
    xlab(expression(paste('Distance to Catalytic Residue (', ring(A), ')'))) +
    geom_smooth(color="black") + 
    geom_vline(xintercept=27.5, size=1, color="red") + 
    coord_cartesian(xlim=c(0,80), ylim=c(0,3))
  
  p1b.label <- (sum(master.clean$dist_active < 27.5)/nrow(master.clean))*100
  
  p1b <- ggplot(master.clean, aes(x=dist_active)) +
    geom_vline(xintercept=shell_breaks, linetype="dashed") +
    geom_density(fill="gray80") + 
    geom_vline(xintercept=27.5, color="red", size=1) +
    xlab(expression(paste('Distance to Catalytic Residue (', ring(A), ')'))) +
    ylab(expression(paste('Density (1/', ring(A), ')'))) + 
    annotate("text", label = paste0(sprintf("%.0f", p1b.label), "%"), x = 18, y = 0.01, size = 7) +
    coord_cartesian(xlim=c(0,80))
  
  p1c <- ggplot(master.clean, aes(x=factor(shell), y=rate)) + 
    geom_violin(scale="width", trim=TRUE, aes(fill=..count..)) +
    stat_summary(fun.y=mean, geom="point") +
    stat_summary(fun.y=mean, geom="line", aes(group=1)) +
    geom_vline(xintercept=c(6.5), color="red", size=1) +
    coord_cartesian(ylim=c(0,3)) +
    scale_x_discrete(limits=as.character(0:10)) +
    xlab('Shell') +
    ylab('Relative Rate') +
    scale_fill_gradient(name="Residue Count", low="lightblue", high="blue", trans="sqrt") +
    theme(legend.key.height = unit(25, "pt"))
  
  # Compile Figure 1 
  fig1 <- ggdraw() +
    draw_plot(p1a, x=0, y=0.5, width=0.5, height=0.5) +
    draw_plot(p1b, x=0.5, y=0.5, width=0.5, height=0.5) +
    draw_plot(p1c, x=0, y=0, width=0.651, height=0.5) +
    draw_plot_label(c('a','b','c'), x = c(0, 0.5, 0), y= c(1,1,0.5))
  
  fig1
  
  # =============================================================
  # FIGURE 2: WCN and RSA
  # =============================================================
  
  # Define color pallete for residual and prediction plots
  cbpallete <- c("alt"= "#E69F00", "d" = "#CC79A7", "null"= "#009E73", "rate" = "black")
  shapes <- c("alt"= 17, "d" = 15, "null"= 16, "rate" = 21)
  # Some preprocessing first...
  
  # Gather predictions for plotting
  master.preds <- select(master.clean, pdb, rate, shell, ends_with('.pred')) %>%
    gather(pred.type, pred, -pdb, -shell) %>%
    group_by(shell, pred.type) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model = get_null(pred.type))
  
  # Gather residuals
  master.res <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  # Fig 2A
  
  # Square all of the correlations and make tidy
  master.lmcor.tidy <- mutate_each(master.lmcor, funs(sq = (.)^2), ends_with('.pred')) %>%
    gather(input, cor, -pdb, -multimer)
  p2a.labels <- c('wcnSC.pred'=expression(K%~%WCN),
                 'RSA.pred'=expression(K%~%RSA),
                 'dist_active.pred'=expression(K%~%d),
                 'wcnSC.RSA.pred'=expression(K%~%WCN+RSA),
                 'wcnSC.dist_active.pred'=expression(K%~%WCN+d), 
                 'RSA.dist_active.pred'=expression(K%~%RSA+d),
                 'wcnSC.RSA.dist_active.pred'=expression(K%~%WCN+RSA+d))
  p2a <- ggplot(master.lmcor.tidy, aes(x=factor(input, c('dist_active.pred',
                                                        'RSA.pred',
                                                        'wcnSC.pred',
                                                        'wcnSC.RSA.pred',
                                                        'wcnSC.dist_active.pred',
                                                        'RSA.dist_active.pred', 
                                                        'wcnSC.RSA.dist_active.pred')), y=cor))  +
    geom_violin(fill="grey80") + 
    scale_x_discrete(labels=p2a.labels) +
    stat_summary(aes(group=input, y=cor, label=sprintf('%0.2f', ..y..)), fun.y=mean, geom="point") + 
    stat_summary(aes(group=input, y=cor, label=sprintf('%0.2f', ..y..)), fun.y=mean, geom="text", hjust=-0.3,vjust=0.5) +
    xlab('Model') + 
    ylab(expression(R^"2")) +
    coord_flip()
  
  # Fig 2B
  cbpallete3 <- c("rate" = "black", "wcnSC.RSA.pred" = "#009E73", "wcnSC.RSA.dist_active.pred"="#E69F00")
  shapes3 <- c("rate" = 21, "wcnSC.RSA.pred" = 16, "wcnSC.RSA.dist_active.pred"=17)
  target <- c('rate','wcnSC.RSA.pred', 'wcnSC.RSA.dist_active.pred')
  p2b <- ggplot(filter(master.preds, pred.type %in% target), 
                 aes(x=factor(shell), y=mean_pred, color=pred.type, shape=pred.type, group=pred.type)) +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-0.1,1.6) +
    labs(x="Shell", y="Relative Rate") +
    scale_color_manual(name="Model", breaks=target, labels=c("Empirical Rate", "WCN + RSA", "WCN + RSA + d"), values=cbpallete3) +
    scale_shape_manual(name="Model", breaks=target, labels=c("Empirical Rate", "WCN + RSA", "WCN + RSA + d"), values=shapes3) +
    theme(legend.justification=c(1,0), legend.position=c(1,-0.02))
  
  # Write out raw data for Fig 2B
  data_out <- filter(master.preds, pred.type %in% target) %>% select(shell, pred.type, mean_pred)
  write_csv(data_out, paste0('models', as.character(suffix),'_data.csv'))
  
  # Fig 2C
  cbpallete2 <- c('wcnSC.RSA.dist_active.res' = "#E69F00", 'wcnSC.RSA.res' = "#009E73")
  shapes2 <- c('wcnSC.RSA.dist_active.res' = 17, 'wcnSC.RSA.res' = 16)
  target <- c('wcnSC.RSA.res', 'wcnSC.RSA.dist_active.res')
  p2c <- ggplot(filter(master.res, pred.type %in% target), 
                aes(x=factor(shell), y=mean_pred, group=pred.type, shape=pred.type, color=pred.type)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-0.6,0.2) +
    labs(x="Shell", y="Residual (Relative Rate)") +
    scale_color_manual(name="Model", breaks=target, labels=c('WCN + RSA', 'WCN + RSA + d'), values=cbpallete2) +
    scale_shape_manual(name="Model", breaks=target, labels=c('WCN + RSA', 'WCN + RSA + d'), values=shapes2) +
    theme(legend.justification=c(1,0), legend.position=c(1,-0.02))
  
  # Compile Figure 2
  fig2 <- ggdraw() +
    draw_plot(p2a, 0, 0, 0.5, 1) +
    draw_plot(p2b, 0.508, 0.5, 0.492, 0.5) +
    draw_plot(p2c, 0.5, 0.009, 0.5, 0.5) +
    draw_plot_label(c('a','b','c'), c(0,0.5,0.5), c(1, 1, 0.53))
  
  # =============================================================
  # FIGURE 3: INDIVIDUAL EXAMPLES
  # =============================================================
  
  # Compute true pearson correlations, without using a linear model
  master.cor <- summarize(master.clean, multimer=max(multimer), 
                          dist_wcn=calc_cor(data=data.frame(y=dist_active,x1=wcnSC)),
                          dist_rsa=calc_cor(data=data.frame(y=dist_active,x1=RSA)),
                          rsa_wcn=calc_cor(data=data.frame(y=RSA,x1=wcnSC)),
                          r4s_dist=calc_cor(data=data.frame(y=rate,x1=dist_active)), 
                          r4s_wcn=calc_cor(data=data.frame(y=rate,x1=wcnSC)), 
                          r4s_rsa=calc_cor(data=data.frame(y=rate,x1=RSA)),
                          r4s_wcnca=calc_cor(data=data.frame(y=rate,x1=wcnCA)))
  
  # We must manually select proteins, otherwise it will be a different set depending on how we filter the data
  # master.cor.ex <- mutate_each(master.cor, funs(abs), -pdb, -multimer) %>%
  #  filter(r4s_dist > 0.5, dist_wcn < 0.25, dist_rsa < 0.25) 
  
  master.cor.ex <- mutate_each(master.cor, funs(abs), -pdb, -multimer) %>% 
    filter(pdb == "1DHF" | pdb == "1DO6" | pdb == "1L0O" | pdb == "1SMN")
  master.ex <- master.clean[master.clean$pdb %in% master.cor.ex$pdb,]
  master.cor.ex <- master.cor.ex %>% 
    select(pdb, r4s_dist, r4s_wcn, r4s_rsa, dist_wcn, dist_rsa, rsa_wcn) %>%
    gather(pair, cor, -pdb)
  p3a <- ggplot(master.cor.ex, aes(x=pair, y=cor, fill=pair)) + 
    geom_bar(stat="identity") +
    scale_x_discrete(labels=c('rsa_wcn' = 'RSA, WCN', 
                              'dist_rsa' = 'd, RSA', 
                              'dist_wcn' = 'd, WCN', 
                              'r4s_rsa' = 'K, RSA', 
                              'r4s_wcn' = 'K, WCN', 
                              'r4s_dist' = 'K, d')) +
    geom_text(aes(group=pair, y=cor, label=sprintf("%0.2f", cor)), hjust=-0.1) +
    xlab('Correlation') +
    ylab(expression('|r|')) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), lim=c(0, 1.1)) +
    facet_wrap(~ pdb, ncol = 1) +
    theme(legend.position = 'none', axis.title.x=element_text(vjust=0.2)) +
    coord_flip() + panel_border()
  
  p3b <- ggplot(master.ex, aes(x=dist_active, y=gam_rate)) +
    geom_point() +
    xlab(expression(paste('Distance to Catalytic Residue (', ring(A), ')'))) +
    ylab("Relative Rate K") +
    facet_wrap(~ pdb, ncol=1) + panel_border()
  
  fig3 <- plot_grid(p3a, p3b, ncol=2, align = 'h', rel_widths = c(1, 0.8), labels=c('a','b'))
  
  fig3
  
  # =============================================================
  # FIGURE 4: ACTIVE SITE LOCATION
  # =============================================================
  
  master.res_fig4 <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type, rsa.group) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  fig4 <- ggplot(filter(master.res_fig4, pred.type != 'rate', pred.type != "dist_active.res"), 
                aes(x=factor(shell), y=mean_pred, group=pred.type, shape=factor(null.model), color=factor(null.model))) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-1.1,0.3) +
    labs(x="Shell", y="Residual (Relative Rate)") +
    scale_color_manual(name="Model", breaks=c('null', 'alt'), labels=c('Structure', 'Structure + d'), values=cbpallete) +
    scale_shape_manual(name="Model", breaks=c('null', 'alt'), labels=c('Structure', 'Structure + d'), values=shapes) +
    facet_grid(null.group ~ rsa.group, labeller=labeller(null.group = label_null_group, rsa.group = label_rsa_group)) +
    panel_border()
  
  # Write raw data for Fig 4 to file
  data_out <- filter(master.res_fig4, pred.type != 'rate', pred.type != "dist_active.res") %>% 
    select(shell, pred.type, mean_pred) %>%
    rename(mean_res = mean_pred)
  write_csv(data_out, paste0('rsa_bins', as.character(suffix),'_data.csv'))
  
  fig4
  
  # =============================================================
  # FIGURE 5: SIZE WITH RSA, WCN, AND DISTANCE
  # =============================================================
  
  # Size vs distance
  p5a <- ggplot(master.clean, aes(x=factor(shell), y=rate)) +
    geom_violin(scale="width", trim=TRUE, aes(fill=..count..)) +
    stat_summary(fun.y=mean, geom="point") +
    stat_summary(fun.y=mean, geom="line", aes(group=1)) +
    coord_cartesian(ylim=c(0,3)) +
    xlab('Shell') +
    ylab('Relative Rate') +
    xlim("0","1","2","3","4","5","6") +
    facet_grid( . ~ ntile.size, labeller = labeller(ntile.size = label_size_group)) +
    scale_fill_gradient(name="Residue Count", low="lightblue", high="blue", trans="sqrt") +
    panel_border()
  
  master.res3 <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type, ntile.size) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  target <- c('wcnSC.RSA.res', 'wcnSC.RSA.dist_active.res')
  p5b <- ggplot(filter(master.res3, pred.type %in% target), 
                 aes(x=factor(shell), y=mean_pred, color=pred.type, group=pred.type, shape=pred.type)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim("0","1","2","3","4","5","6") +
    ylim(-0.6, 0.3) +
    facet_grid(. ~ ntile.size, labeller=labeller(ntile.size = label_size_group)) +
    xlab('Shell') +
    ylab('Residual (Relative Rate)') +
    scale_color_manual(name="Model", breaks=target, labels=c('WCN + RSA', 'WCN + RSA + d'), values=cbpallete2) +
    scale_shape_manual(name="Model", breaks=target, labels=c('WCN + RSA', 'WCN + RSA + d'), values=shapes2) +
    panel_border()
    
  fig5 <- plot_grid(p5a, p5b, ncol = 1, align='v', labels=c('a','b'))
  
  # Write out raw data for Figure 5
  data_out <- filter(master.res3, pred.type %in% target) %>% 
    select(shell, pred.type, mean_pred) %>%
    rename(mean_res = mean_pred)
  write_csv(data_out, paste0('size_bins', as.character(suffix),'_data.csv'))
    
  fig5
  
  # =============================================================
  # FIGURE 8: RSA vs WCN
  # =============================================================
  
  fig8 <- qplot(abs(master.cor$r4s_rsa), 
               abs(master.cor$r4s_wcn), 
               xlab="Correlation (K, RSA)", 
               ylab="Correlation (K, WCN)",
               xlim=c(0:1),
               ylim=c(0:1),
               alpha=0.5) + 
    geom_abline(intercept=0,slope=1) +
    coord_equal() +
    annotate("text", x=0.1, y=0.9, label=sum(abs(master.cor$r4s_wcn)>abs(master.cor$r4s_rsa)), fontface=2) +
    annotate("text", x=0.9, y=0.1, label=sum(abs(master.cor$r4s_wcn)<abs(master.cor$r4s_rsa)), fontface=2) +
    theme(legend.position="none")
  fig8
  
  # =============================================================
  # FIGURE 9: Correlation distributions
  # =============================================================
  
  master.cor.tidy <- master.cor %>%
    select(pdb, multimer, dist_wcn, dist_rsa, rsa_wcn) %>%
    gather(lm, cor, -pdb, -multimer)
  fig9 <- ggplot(master.cor.tidy, aes(x=cor, fill=lm)) + 
    geom_density(alpha=0.25) + 
    labs(x="r", y="Density (1/r)") + 
    scale_fill_manual(name="Combination", 
                      values=c("red","green","blue"), 
                      labels=c("d, WCN", "d, RSA", "WCN, RSA" )) +
    scale_x_continuous(breaks=seq(-1, 1, by=0.25), limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  fig9
  
  # =============================================================
  # FIGURE 10: Model performance
  # =============================================================
  
  master.pred_fig10 <- select(master.clean, pdb, shell, rsa.group, ntile.size, rate, ends_with('.pred')) %>%
    mutate(rate1 = rate, rate2 = rate, rate3 = rate, dist_active.pred2 = dist_active.pred, dist_active.pred3 = dist_active.pred) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  master.res_fig10 <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    mutate(dist_active.res2 = dist_active.res, dist_active.res3 = dist_active.res) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  # Write out raw data for Figure 10a
  data_out <- filter(master.pred_fig10, pred.type != "rate") %>% 
    select(shell, pred.type, mean_pred) 
  write_csv(data_out, paste0('models_dist_line_a', as.character(suffix),'_data.csv'))
  # Write out raw data for Figure 10b
  data_out <- filter(master.res_fig10, pred.type != "rate") %>% 
    select(shell, pred.type, mean_pred) %>%
    rename(mean_res = mean_pred)
  write_csv(data_out, paste0('models_dist_line_b', as.character(suffix),'_data.csv'))
  
  fig10a <- ggplot(filter(master.pred_fig10, pred.type != "rate"), 
                aes(x=factor(shell), y=mean_pred, color=null.model, group=pred.type, shape=null.model)) +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-0.2,2) +
    labs(x="Shell", y="Relative Rate") +
    scale_color_manual(name="Model", breaks=c("rate", "d", "null", "alt"), labels=c("Empirical Rate", "d", "Structure", "Structure + d"), values=cbpallete) +
    scale_shape_manual(name="Model", breaks=c("rate", "d", "null", "alt"), labels=c("Empirical Rate", "d", "Structure", "Structure + d"), values=shapes) +
    facet_grid(null.group ~ ., labeller=labeller(null.group = label_null_group)) + panel_border()
  
  fig10b <- ggplot(filter(master.res_fig10, pred.type != 'rate'), 
                  aes(x=factor(shell), y=mean_pred, group=pred.type, color=factor(null.model), shape=factor(null.model))) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-0.8,0.3) +
    labs(x="Shell", y="Residual (Relative Rate)") +
    scale_color_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d', 'Structure', 'Structure + d'), values=cbpallete) +
    scale_shape_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d', 'Structure', 'Structure + d'), values=shapes) +
    facet_grid(null.group ~ ., labeller=labeller(null.group = label_null_group)) + panel_border()
  
  fig10 <- plot_grid(fig10a, fig10b, labels = c("a","b"))
  
  # =============================================================
  # FIGURE 11: Size grid
  # =============================================================
  
  master.res_fig11 <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    mutate(dist_active.res2 = dist_active.res, dist_active.res3 = dist_active.res) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type, ntile.size) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  fig11 <- ggplot(filter(master.res_fig11, pred.type != 'rate'), 
                 aes(x=factor(shell), y=mean_pred, group=pred.type, color=factor(null.model), shape=factor(null.model))) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-0.7,0.4) +
    labs(x="Shell", y="Residual (Relative Rate)") +
    scale_color_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d', 'Structure', 'Structure + d'), values=cbpallete) +
    scale_shape_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d', 'Structure', 'Structure + d'), values=shapes) +
    facet_grid(null.group ~ ntile.size, labeller=labeller(null.group = label_null_group, ntile.size = label_size_group)) +
    panel_border()
  
  # Write out raw data for Figure 11
  data_out <- filter(master.res_fig11, pred.type != "rate") %>% 
    select(shell, pred.type, mean_pred, ntile.size) %>%
    rename(mean_res = mean_pred, size_group = ntile.size)
  write_csv(data_out, paste0('size_bins_dist_line', as.character(suffix),'_data.csv'))
  
  fig11
  
  # =============================================================
  # FIGURE 12: ACTIVE SITE LOCATION
  # =============================================================
  
  master.res_fig12 <- select(master.clean, pdb, shell, rsa.group, ntile.size, ends_with('.res')) %>%
    mutate(dist_active.res2 = dist_active.res, dist_active.res3 = dist_active.res) %>%
    gather(pred.type, pred, -ntile.size, -pdb, -shell, -rsa.group) %>%
    group_by(shell, pred.type, rsa.group) %>%
    summarize(mean_pred = mean(pred)) %>%
    mutate(null.group = get_null_group(pred.type), null.model=get_null(pred.type))
  
  fig12 <- ggplot(filter(master.res_fig12, pred.type != 'rate'), 
                 aes(x=factor(shell), y=mean_pred, group=pred.type, color=factor(null.model), shape=factor(null.model))) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point() +
    geom_line() +
    xlim(as.character(0:6)) +
    ylim(-1.1,0.3) +
    labs(x="Shell", y="Residual (Relative Rate)") +
    scale_color_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d','Structure', 'Structure + d'), values=cbpallete) +
    scale_shape_manual(name="Model", breaks=c('d', 'null', 'alt'), labels=c('d', 'Structure', 'Structure + d'), values=shapes) +
    facet_grid(null.group ~ rsa.group, labeller=labeller(null.group = label_null_group, rsa.group = label_rsa_group)) +
    panel_border()
  
  # Write out raw data for Figure 11
  data_out <- filter(master.res_fig12, pred.type != "rate") %>% 
    select(shell, pred.type, mean_pred, rsa.group) %>%
    rename(mean_res = mean_pred)
  write_csv(data_out, paste0('rsa_bins_dist_line', as.character(suffix),'_data.csv'))
  
  fig12
  
  # =============================================================
  # SAVE FIGURES
  # =============================================================
  
  save_plot(paste0('distance_rate', as.character(suffix),'.pdf'), fig1, base_width = 10, base_height = 8)
  save_plot(paste0('models', as.character(suffix),'.pdf'), fig2, base_width = 12, base_height = 6)
  save_plot(paste0('examples', as.character(suffix),'.pdf'), fig3, ncol=2, nrow=1, base_height = 10, base_width = 4)
  save_plot(paste0('rsa_bins', as.character(suffix),'.pdf'), fig4, base_width=8, base_height = 6.5)
  save_plot(paste0('size_bins', as.character(suffix),'.pdf'), fig5, ncol=1, base_width = 10, base_height = 6.5)
  save_plot(paste0('rsa_v_wcn', as.character(suffix),'.pdf'), fig8)
  save_plot(paste0('cor_dist', as.character(suffix),'.pdf'), fig9, base_width = 6, base_height = 4)
  save_plot(paste0('models_dist_line', as.character(suffix),'.pdf'), fig10, base_width = 10, base_height = 6.5)
  save_plot(paste0('size_bins_dist_line', as.character(suffix),'.pdf'), fig11, base_width = 8, base_height = 6.5)
  save_plot(paste0('rsa_bins_dist_line', as.character(suffix),'.pdf'), fig12, base_width = 8, base_height = 6.5)
}

# Read in master data table
master <- read_csv("../master_data_table.csv.gz", col_names=TRUE)

# Remove proteins without active site information
master.clean <- group_by(master, pdb,Chain) %>% filter(!is.na(dist_active))

# Define shell breaks
shell_breaks <- c(0, seq(2.5, max(master.clean$dist_active), by=5), max(master.clean$dist_active))

master.clean <- mutate(master.clean, length = n()) %>% # Add length of protein
  ungroup() %>%
  # Cut residues into shells as defined by shell_breaks above
  mutate(shell = cut(dist_active, labels=F, breaks = shell_breaks, include.lowest=T)) %>%
  mutate(shell = shell - 1) %>%
  group_by(pdb)

# Create quartiles for splitting data set based on protein size
quarts <- summarize(master.clean, length=unique(length)) %>% 
  mutate(ntile.size=ntile(length, 3)) %>%
  select(pdb, ntile.size)
# Append quartile data
master.clean <- inner_join(master.clean, quarts, by = c('pdb'='pdb'))

# Compute average RSA for active sites
master.act <- filter(master.clean, ACTIVE_SITE == 1) %>%
  summarize(active_site_rsa = mean(RSA)) %>%
  select(pdb, active_site_rsa)

# Define RSA breaks
# rsa_breaks <- quantile(master.act$active_site_rsa, probs=c(0, 0.333, 0.677, 1))
rsa_breaks <- c(0, 0.05, 0.25, 1)
master.act <- mutate(master.act,
                     rsa.group=cut(active_site_rsa, 
                                   breaks=rsa_breaks, 
                                   include.lowest = T, 
                                   labels=seq(1,3)))
master.clean <- inner_join(master.clean, select(master.act, pdb, rsa.group), by = c('pdb'='pdb'))

# Initial figures with all residues and subunits included
make_figures(master.clean, suffix = "")

# Single subunits, interface included
make_figures(master.clean, MULTI_RSA = F, MULTI_WCN = F, suffix="_mono")

# Single subunits, interface removed
make_figures(master.clean, MULTI_RSA = F, MULTI_WCN = F, NO_INTERFACE = T, suffix="_mono_noint")

# Bio assembly, interface removed
make_figures(master.clean, NO_INTERFACE = T, suffix = "_noint")



# =============================================================
# FIGURE 6 and 7: PREDICTING ACTIVE SITES
# =============================================================

master.act2 <- filter(master.clean, ACTIVE_SITE==1) %>%
  summarize(ref_to_act = min(distance)) %>%
  mutate(group = (ref_to_act < 7.5))
fig6 <- ggplot(master.act2, aes(x=ref_to_act, fill=group)) + 
  geom_histogram(color="black", binwidth=1.5, aes(y=..count../sum(..count..)*100)) +
  scale_fill_manual(values = c("grey80", "red")) +
  xlab(expression(paste("Distance to Catalytic Residue (", ring(A), ")"))) +
  ylab("Enzymes (% of total)") +
  scale_x_continuous(breaks=seq(0,65,10)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 22)) +
  annotate("errorbarh", xmin = 0, xmax = 7.5, y=19, x=0, color="red", size=1, height=1) +
  annotate("text", label=paste0(round(sum(master.act2$ref_to_act<7.5)/nrow(master.act2)*100), "%"), x=3.75, y=20, hjust=0.5, vjust=0) +
  theme(legend.position = "none")
# Write out raw data to produce Figure 6
data_out <- master.act2 %>% select(pdb, group, ref_to_act) %>% rename(near_active = group, ref_to_active = ref_to_act)
write_csv(data_out, "fig_act_data.csv")
fig6

# For WCN control...
master.act4 <- filter(master.clean, wcnSC == max(wcnSC)) %>%
  mutate(group = (dist_active < 7.5))
fig7 <- ggplot(master.act4, aes(x=dist_active, fill=group)) + 
  geom_histogram(color="black", binwidth=1.5, aes(y=..count../sum(..count..)*100)) + 
  scale_fill_manual(values = c("grey80", "red")) +
  xlab(expression(paste("Distance to Catalytic Residue (", ring(A), ")"))) +
  ylab("Enzymes (% of total)") +
  scale_x_continuous(breaks=seq(0,65,10)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 20)) +
  annotate("errorbarh", xmin = 0, xmax = 7.5, y=15, x=0, color="red", size=1, height=1) +
  annotate("text", label=paste0(round(sum(master.act4$dist_active<7.5)/nrow(master.act4)*100), "%"), x=3.75, y=16, hjust=0.5, vjust=0) +
  theme(legend.position = "none")
# Write out raw data to produce Figure 7
data_out <- master.act4 %>% select(pdb, group, dist_active) %>% rename(near_active = group, ref_to_active = dist_active)
write_csv(data_out, "fig_act_control_data.csv")
fig7

# Fisher test for significance
gt7.method2 <- sum(master.act2$ref_to_act >= 7.5)
lt7.method2 <- sum(master.act2$ref_to_act < 7.5)
gt7.method4 <- sum(master.act4$dist_active >= 7.5)
lt7.method4 <- sum(master.act4$dist_active < 7.5)
counts <- matrix(c(lt7.method2, gt7.method2, lt7.method4, gt7.method4), nrow = 2, dimnames = list(c("<7.5", ">7.5"), c("K ~ d", "Max WCN")))
counts
fisher.test(counts)

# Size groups
master.clean %>% 
  group_by(ntile.size) %>% 
  summarize(min_size = min(length), max_size = max(length), count = length(unique(pdb)))

# RSA Groups
master.clean %>% group_by(rsa.group) %>% summarize(count = length(unique(pdb)))

# Monomers/multimers
master.clean %>% ungroup() %>% group_by(multimer) %>% summarize(count = length(unique(pdb)))

save_plot('fig_act.pdf', fig6) # Active site pred
save_plot('fig_act_control.pdf', fig7) # Active site pred with WCN

# system("~/Box\\ Sync/scratch/embed_fonts.py *.pdf")

