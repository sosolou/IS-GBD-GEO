rm(list = ls())
library(easyGBDR)
library(dplyr)
GBD_edition(2021)
library(tidyr)
library(ggplot2)
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/2.overview')

?GBDread
data <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/all.rds") %>%  
  filter(metric!="Percent")


rt <- data %>% 
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>%
  filter(year==2021) %>%
  filter(metric != "Percent") %>% 
  mutate(val_95UI = case_when(
    metric=='Rate'~paste0(sprintf("%.2f", round(val,2)),'\n','(',paste(sprintf("%.2f", round(lower,2)),sprintf("%.2f", round(upper,2)),sep=","),')'),
    metric=='Number'~paste0(sprintf("%.2f", round(val/10000,2)),'\n','(',paste(sprintf("%.2f", round(lower/10000,2)),sprintf("%.2f", round(upper/10000,2)),sep=","),')'))) %>%
  mutate(measure = if_else(measure == 'DALYs (Disability-Adjusted Life Years)', 'DALYs', measure)) %>%
  mutate(measure_sex = paste(measure,sex,sep='_')) %>%
  dplyr::select(location,age,metric,measure_sex,val_95UI) %>%
  filter(!(age=='All ages' & metric=='Rate')) %>%
  filter(!(location!='Global' & (!age %in% c('Age-standardized','All ages')))) %>%
  pivot_wider(names_from = measure_sex, values_from = val_95UI) 
unique(rt$age)


EC <- data %>% 
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>%
  filter(val>0) %>% 
  arrange(measure,location,age,sex,year)  



EC_ASR <- EC %>% filter(age=="Age-standardized")



 ASR_results <- GBDASR_aapc(
  data=EC_ASR,
  startyear = 1990,
  endyear = 2021,
  model = "ln",
  joinpoints = 4,
  rei_included = F,
  CI = TRUE,
  digits = 2,
  sep = " to ",
  constant_variance = F,
  AAPCrange = NULL
)
saveRDS(ASR_results,"ASR_results.rds")


ASR_aapc <- ASR_results$AAPC %>%
  mutate(measure = if_else(measure == 'DALYs (Disability-Adjusted Life Years)', 'DALYs', measure)) %>%
  mutate(measure_sex = paste(measure,sex,sep='_')) %>%
  mutate(AAPC_95CI = case_when(
    P.Value < 0.05  ~ paste0(AAPC_95CI, "*"),
    TRUE ~ AAPC_95CI
  )) %>%
  dplyr::select(location,age,measure_sex,AAPC_95CI) %>%
  pivot_wider(names_from = measure_sex, values_from = AAPC_95CI) %>%
  mutate(metric="AAPC")


order1 <- c("All ages", "Age-standardized")
order2 <- c("Global", "High SDI", "High-middle SDI", "Middle SDI","Low-middle SDI","Low SDI")
order3 <- c("Number","Rate","AAPC")



library(reshape2)
library(ggrepel)
library(tidyverse)
library(dplyr)
region <- data %>%
  filter(!location %in% c("High SDI", "High-middle SDI", 
                          "Middle SDI", "Low-middle SDI",
                          "Low SDI","Global")) %>%
  filter( metric=='Rate' & age == 'Age-standardized') %>%
  filter(measure %in% c('Incidence','Deaths'))  
length(unique(region$location))

data(SDI2021)
force(SDI2021)

SDI <- SDI2021[SDI2021$year >= 1990 & SDI2021$year <= 2021, ] %>% 
  filter(location %in% unique(region$location))
length(unique(SDI$location))

SDI_not_region <- SDI %>% filter(!location %in% unique(region$location))
region_not_SDI <- region %>% filter(!location %in% unique(SDI$location))


order <- c("Low SDI","Low-middle SDI","Middle SDI","High-middle SDI","High SDI")


region_SDI <- region %>%  left_join(SDI,by=c("location","year")) %>% 
  mutate(SDI_class=case_when(
    SDI>=0        & SDI<0.454743 ~ 'Low SDI',
    SDI>=0.454743 & SDI<0.607679 ~ 'Low-middle SDI',
    SDI>=0.607679 & SDI<0.689504 ~ 'Middle SDI',
    SDI>=0.689504 & SDI<0.805129 ~ 'High-middle SDI',
    SDI>=0.805129 & SDI<=1 ~ 'High SDI',
    TRUE ~ ""
  )) %>% 
  mutate(SDI_class = factor(SDI_class, levels = order)) %>% 
  dplyr::select(measure,location,year,val,SDI,SDI_class,everything()) %>% 
  filter(year==2021) 

data3 <- region_SDI %>% filter(measure=="Incidence",sex=="Both")

sdi3 <- ggplot(data = data3, aes(SDI, val, label = location)) +
  geom_point(aes(color = SDI_class)) +
  geom_text_repel(aes(color = SDI_class), size = 2, fontface = 'bold', max.overlaps = 160) +
  geom_smooth(colour='black', stat = "smooth", method = 'loess', se = F, span = 0.5) +
  labs(y = "Incidence in both (ASR, per 100,000)", color = "") +
  labs(x = paste0("SDI (",2021,")"))+
  theme(legend.position = c(0.09, 0.80), legend.text = element_text(size = 7), 
        legend.background = element_blank(), axis.text = element_text(size = 7),
        axis.title = element_text(size = 9)) +
  annotate("text", x = min(data3$SDI), y = max(data3$val), 
           label = paste0("Spearman r =", round(cor(data3$SDI, data3$val, method="spearman"), 2), ", P ", ifelse(cor.test(data3$SDI, data3$val, method="spearman")$p.value < 0.001, "<0.001", 
                                                                                                                 ifelse(cor.test(data3$SDI, data3$val, method="spearman")$p.value < 0.01, "<0.01", 
                                                                                                                        ifelse(cor.test(data3$SDI, data3$val, method="spearman")$p.value < 0.05, "<0.05", 
                                                                                                                               ">0.5")))), 
           size = 3, hjust = 0, vjust = 1)

library(gridExtra)
png("region_sdi_ASR_2021.png", width = 10, height = 6, units = "in", res = 300)
grid.arrange(sdi3,ncol = 1)
dev.off()



top10_country.incidence <- data3 %>% top_n(10, val) %>% arrange(desc(val))


top10_country.death <- data3 %>% top_n(10, val) %>% arrange(desc(val))


order <- c("Global", "High SDI", "High-middle SDI", "Middle SDI","Low-middle SDI","Low SDI")

ASR_apc <- ASR_results$APC


ASR_aapc <- ASR_results$AAPC %>% 
  mutate(location = factor(location, levels = order)) 

global_aapc <- ASR_aapc %>% 
  filter(location=="Global") %>% 
  filter(measure %in% c("Incidence","Deaths"))

ASR_aapc1 <- ASR_aapc %>% 
  mutate(AAPC_95CI = case_when(
    P.Value < 0.05  ~ gsub("\n", "",paste0(AAPC_95CI, "*")),
    TRUE ~ gsub("\n", "",AAPC_95CI) 
  )) %>% 
  dplyr::select(measure,location,sex,AAPC_95CI,Test.Statistic,P.Value,everything()) %>% 
  arrange(measure,sex,location)


unique(ASR_aapc1$location)

#######不同SDI地区、不同性别#######
Incidence_global <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASIR",
  point_color = "black",
  joinpoint_color = "black") + 
  #scale_y_continuous(limits = c(0,3)) + #设定y轴的最大值和最小值
  labs(y="ASIR (per 100,000)")  +#修改X轴标签
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))

Incidence_male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASIR for Male",
  point_color = "black",
  joinpoint_color = "black") + 
  #scale_y_continuous(limits = c(0,3)) + #设定y轴的最大值和最小值
  labs(y="ASIR (per 100,000)")  +#修改X轴标签
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))

Incidence_female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASIR for Female",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))

Deaths_global <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Deaths",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASMR",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASMR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))

Deaths_male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Deaths",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASMR for Male",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASMR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))

Deaths_female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Global",
  measure_name = "Deaths",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "Global IS ASMR for Female",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASMR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.1,0),
        text = element_text(size = 24),
        axis.text = element_text(size = 24),
        facet.title.text=element_text(size = 24))


library(gridExtra)
png("joinpoint.incidence and death.png", width = 45, height = 25, units = "in", res = 300)
grid.arrange(Incidence_global,Incidence_male,Incidence_female,
             Deaths_global,Deaths_male,Deaths_female,
             ncol = 3)
dev.off() 


Incidence_High <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR in High SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_High_Male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for male in High SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_High_Female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for Female in High SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))



Incidence_Highmiddle <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR in High-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Highmiddle_Male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for male in High-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Highmiddle_Female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "High-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for Female in High-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_middle <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR in Middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_middle_Male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for male in Middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_middle_Female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for Female in Middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Lowmiddle <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR in Low-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Lowmiddle_Male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for male in Low-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Lowmiddle_Female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low-middle SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for Female in Low-middle SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Low <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Both",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR in Low SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))


Incidence_Low_Male <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Male",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for male in Low SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

Incidence_Low_Female <- ggjoinpoint_apc(
  data = ASR_results,
  location_name = "Low SDI",
  measure_name = "Incidence",
  cause_name = "Ischemic stroke",
  sex_name = "Female",
  age_name = "Age-standardized",
  rei_name = NULL,
  facet_name = "IS ASIR for Female in Low SDI regions",
  point_color = "black",
  joinpoint_color = "black") + 
  labs(y="ASIR (per 100,000)")  +
  theme(legend.position = c(1,0),
        legend.justification = c(3.5,0),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        facet.title.text=element_text(size = 20))

library(gridExtra)
png("joinpoint.incidence.sdi.png", width = 40, height = 40, units = "in", res = 400)
grid.arrange(Incidence_High,Incidence_High_Male,Incidence_High_Female,
             Incidence_Highmiddle,Incidence_Highmiddle_Male,Incidence_Highmiddle_Female,
             Incidence_middle,Incidence_middle_Male,Incidence_middle_Female,
             Incidence_Lowmiddle,Incidence_Lowmiddle_Male,Incidence_Lowmiddle_Female,
             Incidence_Low,Incidence_Low_Male,Incidence_Low_Female,
             ncol = 3)
dev.off() 




