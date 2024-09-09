rm(list = ls())
library(easyGBDR)
library(dplyr)
GBD_edition(2021)
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/3.bapc')
ls("package:easyGBDR")

data <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")
str(data)

order <- c("Global","High SDI","High-middle SDI","Middle SDI","Low-middle SDI","Low SDI")
EC <- data %>% 
  filter(location %in% order) %>%
  filter(val>0) %>% 
  filter(metric != "Percent") %>% 
  arrange(measure,location,sex,age,year)    

unique(EC$measure)
unique(EC$sex)
unique(EC$age)
unique(EC$cause)
unique(EC$metric)
length(unique(EC$year))
unique(EC$location)

bapc_COVID <- GBDbapc_COVID(
  data = EC,
  measure_name = c("Incidence", "Deaths" ),
  cause_name = "Ischemic stroke",
  location_name = unique(EC$location),
  rei_name = NULL,
  full_age_adjusted = T,
  rate_lessen = NULL
)

bapc_COVID_All <-bapc_COVID[["all_age_projection"]] %>% 
  filter(year %in% c(2020,2021)) %>% 
  rename(true_val = ture_val) %>% 
  mutate(cz=round(true_val-expected_val,0),index="Numbers") %>% 
  mutate(true_val=round(true_val,0)) %>% 
  mutate(expected_val=round(expected_val,0)) %>% 
  mutate(location = factor(location, levels = order)) %>% 
  arrange(measure,year,sex,location) %>% 
  dplyr::select(index,measure,year,sex,location,true_val,expected_val,cz)

bapc_COVID_ASR <- bapc_COVID[["ASR"]] %>% 
  filter(year %in% c(2020,2021)) %>% 
  rename(true_val = ture_val) %>% 
  mutate(cz=true_val-expected_val,index="ASR") %>% 
  mutate(location = factor(location, levels = order)) %>% 
  arrange(measure,year,sex,location) %>% 
  dplyr::select(index,measure,year,sex,location,true_val,expected_val,cz)

bapc_COVID_out <- rbind(bapc_COVID_All,bapc_COVID_ASR)
write.csv(bapc_COVID_out,"bapc_COVID.csv",row.names = T)


bapc_result <- GBDbapc_prediction(
  data = EC,
  measure_name = c("Incidence", "Deaths" ),
  cause_name = "Ischemic stroke",
  location_name = "Global",
  rei_name = NULL,
  By_sex = T,
  predyear = 2040,
  full_age_adjusted = T,
  rate_lessen = NULL
)
saveRDS(bapc_result,"bapc_result.rds")



a1 <- bapc_result[["all_age_projection"]] %>%
  filter(year==2040) %>%
  mutate(metric="Number",age="All ages")

a2 <- bapc_result[["crude_rate"]] %>%
  filter(year==2040) %>%
  mutate(metric="Crude-rate",age="All ages")

a3 <- bapc_result[["ASR"]] %>%
  filter(year==2040) %>%
  mutate(metric="ASR",age="Age-standardized")


bapc <- rbind(a1,a2,a3) %>%
  mutate(pred = paste0(round(pred_val,2),' (',
           paste(round(pred_low,2),round(pred_up,2),sep = ' to '),')')
         ) %>%
  dplyr::select(measure,location,age,metric,sex,pred,cause,year) %>%
  arrange(measure,location,metric,age,sex) 

write.csv(bapc,"bapc.csv",row.names = F)



library(ggplot2)

a1 <- ggprediction_Dx(
  bapc_result,
  ratio='auto',
  CI = F,
  predict_start = 2022,
  group_name = "measure",
  location_name = "Global",
  measure_name = c("Incidence","Deaths"),
  cause_name = "Ischemic stroke",
  rei_name = NULL,
  sex_name = c("Both",'Male','Female')
)+ facet_wrap(.~sex)


library(gridExtra)
png("bapc.png", width = 12, height = 5, units = "in", res = 300)
grid.arrange(a1,ncol = 1)
dev.off()




