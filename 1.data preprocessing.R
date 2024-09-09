rm(list = ls()) 
library(easyGBDR)
library(dplyr)
GBD_edition(2021)
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD')



data <- GBDread(foldername = './1.data', 
                folder = T)
str(data)

data1 <- GBDtrans_normalization(data)        
a <- data1[-which(data1$location %in% data$location),]

unique(data$measure)
unique(data$location)
unique(data$sex)
unique(data$age)
unique(data$cause)
unique(data$metric)
unique(data$year)

data0 <- data %>% filter(!age %in% c("Age-standardized","<5 years")) %>% 
  arrange(age,year)

colnames(data0)
unique(data0$measure)
unique(data0$location)
unique(data0$sex)
unique(data0$age)
b <- length(data0$val)
a <- length(unique(data0$measure))*
  length(unique(data0$location))*
  length(unique(data0$sex))*
  length(unique(data0$age))*
  length(unique(data0$cause))*
  length(unique(data0$metric))
b/a


all <- data %>% filter(age %in% c("All ages","Age-standardized"))
all <- all %>% mutate(location=case_when(location=="T眉rkiye"~"Turkey",
                                           TRUE~location))
 
saveRDS(all,"all.rds")
saveRDS(data,"data.rds")


