rm(list = ls())
library(easyGBDR)
GBD_edition(2021)
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/9.map')

data <- readRDS("C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/all.rds") %>% 
  filter(metric!="Percent")
str(data)

unique(data$measure)
unique(data$sex)
unique(data$age)
unique(data$cause)
unique(data$metric)
unique(data$year)
length(unique(data$location))


data("GBDRegion2021")


############1.Incidence map####################
#遍历32年Incidence ASR
years <- 2021
  map_Incidence <- data %>% 
    filter(measure == 'Incidence' & sex == 'Both' & year == years
           & metric == 'Rate' & age == 'Age-standardized') %>%
    filter(location %in% GBDRegion2021$location) %>%
    mutate(val2 = Quantile(val, n=10))
  
  ggGBDmap(
    data = map_Incidence,
    variable = 'val',
    color = paste0("scale_fill_distiller(palette='Spectral',name=paste0('Incidence in ', years, '\n(ASR, per 100,000)')
                   )"),
    guide_name = 'ASR'
  )
  
  ggsave(paste0('Incidence', years, '.png'), width = 10, height = 6.5, dpi = 300)
