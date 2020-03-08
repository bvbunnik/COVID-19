library(tidyverse)

#startday
var_start_day_scen1 <- read_csv("/media/bram/DATA/covid-19/output/scenario1_vary_startday.csv")
min_peak_I_start_day_scen1 <- var_start_day_scen1 %>% 
  group_by(start_day) %>% 
  summarise(maximum=max(I))

ggplot(min_peak_I_start_day_scen1, aes(x=start_day, y=maximum)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Trigger point (days)", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.16),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))



var_start_day_scen3 <- read_csv("/media/bram/DATA/covid-19/output/scenario3_vary_startday.csv")
min_peak_I_start_day_scen3 <- var_start_day_scen3 %>% 
  group_by(start_day) %>% 
  summarise(maximum=max(I))

ggplot(min_peak_I_start_day_scen3, aes(x=start_day, y=maximum)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Trigger point (days)", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.16),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

comb_data_start_day <- rbind(min_peak_I_start_day_scen1, min_peak_I_start_day_scen3)
comb_data_start_day$scen <- rep(c(1,3),times=c(100,100))

ggplot(comb_data_start_day, aes(x=start_day, y=maximum, colour = as.factor(scen))) + geom_line(size = 1.05) +
  labs(x ="Trigger point (days)", y = "Peak fraction infected", colour="Scenario") +
  scale_y_continuous(limits = c(0,0.16),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  
  scale_x_continuous(expand = c(0, 0))


write_csv(comb_data_start_day, "/media/bram/DATA/covid-19/output/Uncertainty_triggerpoint.csv")



#beta
var_beta_scen1 <- read_csv("/media/bram/DATA/covid-19/output/scenario1_varying_beta.csv")
min_peak_I_beta_scen1 <- var_beta_scen1 %>% 
  group_by(beta_scen) %>% 
  summarise(maximum=max(I)) 

ggplot(min_peak_I_beta_scen1, aes(x=beta_scen, y=maximum))+geom_line(size=1.05) +
  labs(x ="Fractional reduction in beta", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.15),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

min_peak_I_beta_scen1$scen <- 1

var_beta_scen3 <- read_csv("/media/bram/DATA/covid-19/output/scenario3_varying_beta.csv")
min_peak_I_beta_scen3 <- var_beta_scen3 %>% 
  group_by(beta_scen) %>% 
  summarise(maximum=max(I)) 


ggplot(min_peak_I_beta_scen3, aes(x=beta_scen, y=maximum))+geom_line(size=1.05) +
  labs(x ="Fractional reduction in beta", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.15),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

min_peak_I_beta_scen3$scen <- 3

comb_data_beta <- rbind(min_peak_I_beta_scen1, min_peak_I_beta_scen3)
comb_data_beta$scen <- as_factor(comb_data1$scen)

ggplot(comb_data_beta, aes(x=beta_scen, y=maximum, colour=scen))+geom_line(size=1.05) +
  labs(x ="Fractional reduction in beta", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.15),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_discrete(name="Scenario")

write_csv(comb_data_beta, "/media/bram/DATA/covid-19/output/Uncertainty_beta.csv")


#R0
var_R0_scen1 <- read_csv("/media/bram/DATA/covid-19/output/scenario1_varying_R0_startday.csv")
min_peak_I_R0_scen1 <- var_R0_scen1 %>% 
  group_by(R0) %>% 
  summarise(maximum=max(I)) 

ggplot(min_peak_I_R0_scen1, aes(x=R0, y=maximum))+geom_line(size=1.05) +
  labs(x ="R0", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.57),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

min_peak_I_R0_scen1$scen <- 1

var_R0_scen3 <- read_csv("/media/bram/DATA/covid-19/output/scenario3_varying_R0_startday.csv")
min_peak_I_R0_scen3 <- var_R0_scen3 %>% 
  group_by(R0) %>% 
  summarise(maximum=max(I)) 


ggplot(min_peak_I_R0_scen3, aes(x=R0, y=maximum))+geom_line(size=1.05) +
  labs(x ="R0", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.55),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

min_peak_I_R0_scen3$scen <- 3

comb_data_R0 <- rbind(min_peak_I_R0_scen1, min_peak_I_R0_scen3)
comb_data_R0$scen <- as_factor(comb_data1$scen)

ggplot(comb_data_R0, aes(x=R0, y=maximum, colour=scen))+geom_line(size=1.05) +
  labs(x ="R0", y = "Peak fraction infected") +
  scale_y_continuous(limits = c(0,0.31),  expand = c(0,0)) +
  theme(legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_discrete(name="Scenario")

write_csv(comb_data_R0, "/media/bram/DATA/covid-19/output/Uncertainty_R0.csv")