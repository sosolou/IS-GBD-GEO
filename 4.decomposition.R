rm(list = ls())
library(easyGBDR)
GBD_edition(2021)
library(tidyverse)


setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/4.decomposition')

order1 <- c("Global", "High SDI", "High-middle SDI", "Middle SDI","Low-middle SDI","Low SDI")
order2 <- c("Both","Male","Female")

data <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>%
  filter(location %in% order1) %>% 
  mutate(location = factor(location, levels = order1)) %>%
  mutate(sex = factor(sex, levels = order2)) %>%
  arrange(desc(measure),location,sex,age,year)


data1990_1999 <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>%
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>%  
  filter(year %in% c(1990,1999) & age=='All ages' & metric=='Number') %>% 
  dplyr::select(location,sex,measure,year,val) %>%
  pivot_wider(names_from = year, values_from = val) 
 

unique(data$age)

results_1990 <- GBDdecomposition(
  data,
  byear = 1990,
  compareyear = 1999,
  startage = 0,
  endage = 95
)


table_incidence1990 <- GBDdecomposition_table(
  data = results_1990,
  digits = 2,
  measure_name = unique(results_1990$measure),
  location_name = unique(results_1990$location),
  sex_name = unique(results_1990$sex),
  rei_name = NULL,
  cause_name = unique(results_1990$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>%
  filter(measure =="Incidence") %>% 
  left_join(data1990_1999,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"1990","1999",everything()) %>% 
  mutate(percent=paste0(round(`1999`/`1990`*100-100,2),"%"))
write.csv(table_incidence1990,"table_incidence1990.csv",row.names = F)


options(scipen=200) 
incidence_1990 <- ggdecomposition(
  data = results_1990,
  measure_name = "Incidence",
  location_name = unique(results_1990$location),
  sex_name = unique(results_1990$sex),
  rei_name = NULL,
  cause_name = unique(results_1990$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Incidence number changes (1990-1999)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(fill=FALSE)


results_death1990 <- GBDdecomposition_servity(
  data,
  byear = 1990,
  compareyear = 1999,
  measure_name = "Deaths",
  startage = 0,
  endage = 95
)

table_death1990 <- GBDdecomposition_table(
  data = results_death1990,
  digits = 2,
  measure_name = unique(results_death1990$measure),
  location_name = unique(results_death1990$location),
  sex_name = unique(results_death1990$sex),
  rei_name = NULL,
  cause_name = unique(results_death1990$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>% 
  left_join(data1990_1999,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"1990","1999",everything())%>% 
  mutate(percent=paste0(round(`1999`/`1990`*100-100,2),"%"))
write.csv(table_death1990,"table_death1990.csv",row.names = F)

death_1990 <- ggdecomposition(
  data = results_death1990,
  measure_name = "Deaths",
  location_name = unique(results_death1990$location),
  sex_name = unique(results_death1990$sex),
  rei_name = NULL,
  cause_name = unique(results_death1990$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Death number changes (1990-1999)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(fill=FALSE)


data2000_2009 <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>%
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>%  
  filter(year %in% c(2000,2009) & age=='All ages' & metric=='Number') %>% 
  dplyr::select(location,sex,measure,year,val) %>%
  pivot_wider(names_from = year, values_from = val) 


results_2000 <- GBDdecomposition(
  data,
  byear = 2000,
  compareyear = 2009,
  startage = 0,
  endage = 95
)


table_incidence2000 <- GBDdecomposition_table(
  data = results_2000,
  digits = 2,
  measure_name = unique(results_2000$measure),
  location_name = unique(results_2000$location),
  sex_name = unique(results_2000$sex),
  rei_name = NULL,
  cause_name = unique(results_2000$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>%
  filter(measure =="Incidence") %>% 
  left_join(data2000_2009,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"2000","2009",everything()) %>% 
  mutate(percent=paste0(round(`2009`/`2000`*100-100,2),"%"))
write.csv(table_incidence2000,"table_incidence2000.csv",row.names = F)

#?ggdecomposition
options(scipen=200) #不要使用科学计数法
incidence_2000 <- ggdecomposition(
  data = results_2000,
  measure_name = "Incidence",
  location_name = unique(results_2000$location),
  sex_name = unique(results_2000$sex),
  rei_name = NULL,
  cause_name = unique(results_2000$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Incidence number changes (2000-2009)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(fill=FALSE)


results_death2000 <- GBDdecomposition_servity(
  data,
  byear = 2000,
  compareyear = 2009,
  measure_name = "Deaths",
  startage = 0,
  endage = 95
)

table_death2000 <- GBDdecomposition_table(
  data = results_death2000,
  digits = 2,
  measure_name = unique(results_death2000$measure),
  location_name = unique(results_death2000$location),
  sex_name = unique(results_death2000$sex),
  rei_name = NULL,
  cause_name = unique(results_death2000$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>% 
  left_join(data2000_2009,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"2000","2009",everything()) %>% 
  mutate(percent=paste0(round(`2009`/`2000`*100-100,2),"%"))
write.csv(table_death2000,"table_death2000.csv",row.names = F)

death_2000 <- ggdecomposition(
  data = results_death2000,
  measure_name = "Deaths",
  location_name = unique(results_death2000$location),
  sex_name = unique(results_death2000$sex),
  rei_name = NULL,
  cause_name = unique(results_death2000$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Death number changes (2000-2009)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  guides(fill=FALSE)


data2010_2021 <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/data.rds")  %>%
  filter(location %in% c("Global","High SDI", "High-middle SDI", 
                         "Middle SDI", "Low-middle SDI",
                         "Low SDI")) %>%  
  filter(year %in% c(2010,2021) & age=='All ages' & metric=='Number') %>% 
  dplyr::select(location,sex,measure,year,val) %>%
  pivot_wider(names_from = year, values_from = val) 


results_2010 <- GBDdecomposition(
  data,
  byear = 2010,
  compareyear = 2021,
  startage = 0,
  endage = 95
) 

#?GBDdecomposition_table()
table_incidence2010 <- GBDdecomposition_table(
  data = results_2010,
  digits = 2,
  measure_name = unique(results_2010$measure),
  location_name = unique(results_2010$location),
  sex_name = unique(results_2010$sex),
  rei_name = NULL,
  cause_name = unique(results_2010$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>%
  filter(measure =="Incidence") %>% 
  left_join(data2010_2021,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"2010","2021",everything()) %>% 
  mutate(percent=paste0(round(`2021`/`2010`*100-100,2),"%"))
write.csv(table_incidence2010,"table_incidence2010.csv",row.names = F)


options(scipen=200) 
incidence_2010 <- ggdecomposition(
  data = results_2010,
  measure_name = "Incidence",
  location_name = unique(results_2010$location),
  sex_name = unique(results_2010$sex),
  rei_name = NULL,
  cause_name = unique(results_2010$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Incidence number changes (2010-2021)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) 


results_death2010 <- GBDdecomposition_servity(
  data,
  byear = 2010,
  compareyear = 2021,
  measure_name = "Deaths",
  startage = 0,
  endage = 95
)

table_death2010 <- GBDdecomposition_table(
  data = results_death2010,
  digits = 2,
  measure_name = unique(results_death2010$measure),
  location_name = unique(results_death2010$location),
  sex_name = unique(results_death2010$sex),
  rei_name = NULL,
  cause_name = unique(results_death2010$cause)
) %>% 
  dplyr::select(-cause) %>% 
  arrange(measure,location,sex) %>% 
  left_join(data2010_2021,by=c("measure","location","sex")) %>% 
  dplyr::select(location,sex,measure,"2010","2021",everything()) %>% 
  mutate(percent=paste0(round(`2021`/`2010`*100-100,2),"%"))
write.csv(table_death2010,"table_death2010.csv",row.names = F)

death_2010 <- ggdecomposition(
  data = results_death2010,
  measure_name = "Deaths",
  location_name = unique(results_death2010$location),
  sex_name = unique(results_death2010$sex),
  rei_name = NULL,
  cause_name = unique(results_death2010$cause)
) + 
  facet_wrap(.~sex) +
  labs(fill = "Factor",y="Death number changes (2010-2021)",x="")+
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) 


library(gridExtra)
png("decomposition.png", width = 30, height = 20, units = "in", res = 300)
grid.arrange(incidence_1990,death_1990,incidence_2000,death_2000,incidence_2010,death_2010,
             ncol = 2)
dev.off()



