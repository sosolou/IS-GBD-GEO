rm(list = ls())
library(easyGBDR)
GBD_edition(2021)
library(tidyverse)

setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/5.apc')

EC <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>% 
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                          "Middle SDI", "Low-middle SDI",
                          "Low SDI")) %>% 
  filter(val > 0) 

global <- EC %>% 
  filter(location == "Global", sex=="Both",metric=="Rate",measure=="Incidence",year==2021)


apc_web <- GBDapc_web(
  data = EC,
  startyear = 1992,
  endyear = 2021,
  reference_age = 15,
  reference_year = 1992,
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  sex_name = unique(EC$sex),
  location_name = c("Global","High SDI", "High-middle SDI", 
                    "Middle SDI", "Low-middle SDI",
                    "Low SDI"),
  rei_name = NULL
)

  unique(EC$age)
saveRDS(apc_web,file = "apc_web.rds")


drift.global <- ggdrift(
  data = apc_web,
  CI = F,
  group_name = "sex",
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  sex_name = c('Both','Male','Female'),
  location_name = "Global") + facet_wrap(.~measure) +
  labs(title = "Net and local drift") +
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())

effect1 <- ggLongAge(
  data = apc_web,
  CI = F,
  group_name = "sex",
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  sex_name = c('Both','Male','Female'),
  location_name = "Global") + facet_wrap(.~measure)+
  labs(title = "Age effects") +
  theme(plot.title = element_text(hjust = 0.5))

effect2 <- ggcohortRR(
  data = apc_web,
  CI = F,
  group_name = "sex",
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  sex_name = c('Both','Male','Female'),
  location_name = "Global") + facet_wrap(.~measure) +
  labs(title = "Cohort effects") +
  theme(plot.title = element_text(hjust = 0.5))

effect3 <- ggperiodRR(
  data = apc_web,
  CI = F,
  group_name = "sex",
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  sex_name = c('Both','Male','Female'),
  location_name = "Global") + facet_wrap(.~measure)+
  labs(title = "Period effects") +
  theme(plot.title = element_text(hjust = 0.5))

library(gridExtra)
png("APC.global.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(drift.global,effect1,effect2,effect3, ncol = 2)
dev.off()

