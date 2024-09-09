rm(list = ls()) 
library(easyGBDR)
GBD_edition(2021)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)


setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/8.index')
data <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/all.rds") 


rt <- data %>% 
  filter(!(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI"))) %>% 
  filter(age %in% c("All ages","Age-standardized"),year %in% c(1990,2021)) 

length(unique(rt$location))
?GBDslope_index
res <-  GBDslope_index(rt,
                       all_age_range = NULL,
                       SDI = T,
                       GBDregion = T,
                       SuperGBDregion = T)


slop <- res[["slope"]] %>% 
  filter(measure=='DALYs (Disability-Adjusted Life Years)',age=='All ages') %>% 
  dplyr::select(region,sex,year,rlm,everything())%>% 
  arrange(region,sex,year) %>% 
  filter(region %in% c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI","Low SDI","All included"))
slop_out <- slop %>% dplyr::select(region,measure,sex,year,rlm,rlm_lower,rlm_upper)
write.csv(slop_out,"slop_out.csv",row.names = F)                                                                                                                                                                                                           

DALYs_both_slope <- ggslope_index(
  data = res,
  model = "rlm",
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Both",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label = "China" ,
  population_count = 10e6
) + labs(title = "DALYs slope index for both")

DALYs_female_slope <- ggslope_index(
  data = res,
  model = "rlm",
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Female",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label = "China" ,
  population_count = 10e6
) + labs(title = "DALYs slope index for female")

DALYs_male_slope <- ggslope_index(
  data = res,
  model = "rlm",
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Male",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label = "China" ,
  population_count = 10e6
) + labs(title = "DALYs slope index for male")


config_stata(path='D:/ProgramFiles/Stata17',version=17,stata_type='MP')
stata_package_install('conindex')
stata_package_install('apc')
library(ggbrace)
res1 <-  GBDconcentration_index(
  data = rt,
  all_age_range = NULL,
  SDI = T,
  GBDregion = T,
  SuperGBDregion = T
)
saveRDS(res1,file = "res1.rds")

?ggconcentration_index
DALYs_both_concentration <- ggconcentration_index(
  data = res1,
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Both",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label ="China",
  population_count = 1e+06,
  line_type = "geom_line"
)+ labs(title = "DALYs concentration index for both")

DALYs_female_concentration <- ggconcentration_index(
  data = res1,
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Female",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label ="China",
  population_count = 1e+06,
  line_type = "geom_line"
)+ labs(title = "DALYs concentration index for female")

  
DALYs_male_concentration <- ggconcentration_index(
  data = res1,
  color_name = c("#6699FF", "#990000"),
  group_name =  "year",
  region_name = "All included",
  measure_name = "DALYs (Disability-Adjusted Life Years)",
  sex_name = "Male",
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  age_name = "All ages",
  year_name = c(1990, 2021),
  country_label ="China",
  population_count = 1e+06,
  line_type = "geom_line"
)+ labs(title = "DALYs concentration index for male")  


slope <- res[["slope"]] %>% 
  filter(region=="All included", measure=="DALYs (Disability-Adjusted Life Years)") %>% 
  arrange(sex,age,year) %>% 
  dplyr::select(-lm,-lm_lower,-lm_upper,-ncvTest_p)
  
library(gridExtra)
png("DALYs slop and concentration index.png", width = 20, height = 10, units = "in", res = 300)
grid.arrange(DALYs_both_slope,DALYs_female_slope,DALYs_male_slope,
             DALYs_both_concentration,DALYs_female_concentration,DALYs_male_concentration,
             ncol = 3)
dev.off() 

