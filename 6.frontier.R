rm(list = ls())

library(easyGBDR)
GBD_edition(2021)
library(tidyverse)


setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/6.frontier/Both')

case <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>%
  filter(!location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>% 
  filter(age %in% c("All ages", "Age-standardized")) %>% 
  filter(sex=="Both")
length(unique(case$location))


boostrap_DEA_DALYs <- GBDfrontier2(
  data = case,
  sex_name = "Both",
  cause_name = "Ischemic stroke",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  rei_name = NULL,
  age_name = "Age-standardized",
  boot = 10
)

p1 <- ggfrontier(
  boostrap_DEA_DALYs,
  smooth_span = 0.3,
  type = "single year",
           high_SDI = 0.85,
           low_SDI = 0.5
) +
  labs(y="DALYs (ASR, per 100,000)", color="Trend")


boostrap_DEA_DALYs <- GBDfrontier2(
  data = case,
  sex_name = "Female",
  cause_name = "Ischemic stroke",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  rei_name = NULL,
  age_name = "Age-standardized",
  boot = 10
)


p2 <- ggfrontier(
  boostrap_DEA_DALYs,
  smooth_span = 0.3,
  type = "single year",
  high_SDI = 0.85,
  low_SDI = 0.5
) +
  labs(y="DALYs (ASR, per 100,000)", color="Trend")
saveRDS(p2,"p1.rds")

boostrap_DEA_DALYs <- GBDfrontier2(
  data = case,
  sex_name = "Male",
  cause_name = "Ischemic stroke",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  rei_name = NULL,
  age_name = "Age-standardized",
  boot = 10
)


p3 <- ggfrontier(
  boostrap_DEA_DALYs,
  smooth_span = 0.3,
  type = "single year",
  high_SDI = 0.85,
  low_SDI = 0.5
) +
  labs(y="DALYs (ASR, per 100,000)", color="Trend")
