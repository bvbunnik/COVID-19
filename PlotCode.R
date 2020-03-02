rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("/home/bram/Documents/COVID-19/output/final")

create_plot <- function(data, panels=TRUE, limit_inf = c(0,0.20), limit_beta=c(0,0.3)){
  inf_plot <- ggplot(data, aes(x = t, y = I)) + geom_line(size = 1.05, col = "darkblue") +
    labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
    scale_y_continuous(limits = limit_inf ,  expand = c(0,0)) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
    scale_x_continuous(expand = c(0, 0))
  
  if (panels==TRUE){
    inf_plot <- inf_plot + facet_wrap(~scen, nrow=1)
  }
  
  beta_plot <- ggplot(data, aes(x = t, y = beta)) + geom_line(size = 1.05, col = "darkred") +
    labs(x ="Time (Days)", y = "β") +
    scale_y_continuous(limits = limit_beta ,  expand = c(0,0)) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
    scale_x_continuous(expand = c(0, 0)) 
  if (panels==TRUE){
    beta_plot <- beta_plot + facet_wrap(~scen, nrow=1) + theme(strip.background = element_blank(), strip.text = element_blank())
  }

  comb_plot <- ggarrange(inf_plot + rremove("xlab"), beta_plot, nrow=2, heights = c(1, 0.5))
  return(comb_plot)
}

#########################
# Baseline              #
#########################
baseline_all <- read_csv("output_all_scenarios_baseline.csv")
baseline_all$scen = as.factor(baseline_all$scen)
levels(baseline_all$scen) <- c("Baseline", "S1", "S2", "S3", "S4", "S5")
create_plot(baseline_all)
