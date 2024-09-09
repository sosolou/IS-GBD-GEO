rm(list = ls()) 
library(easyGBDR)
GBD_edition(2021)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)

setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/7.risk') 
order1 <- read.csv('order.csv',header = F)
order1$V1 <- rev(order1$V1) 


data <- GBDread(foldername = '.',
                folder = T) %>% 
  filter(measure=="DALYs (Disability-Adjusted Life Years)")
order <- read.csv('order.csv',header = F)
order$V1 <- rev(order$V1) 

Risk_2021 <- subset(data, year==2021 &
                      metric=='Percent' &
                      sex=='Both' &
                      age=='All ages' )
Risk_2021$val <- round(Risk_2021$val*100,1)
Risk_2021$val2 <- paste0(Risk_2021$val,'%')
Risk_2021$location <- factor(Risk_2021$location, 
                             levels=order$V1, 
                             ordered=TRUE)
Risk_2021_1 <- Risk_2021 %>% filter(location=='Global') %>% arrange(desc(measure),desc(val))

factors  <- Risk_2021_1[1:6,]$rei

Risk_2021$rei <- factor(Risk_2021$rei, 
                        levels= factors, 
                        ordered=TRUE)  
Risk_2021 <- na.omit(Risk_2021) 
Risk_2021[which(Risk_2021$measure=='DALYs (Disability-Adjusted Life Years)'),]$measure <-"DALYs"


color <- c("pink2","lightskyblue3","plum3","salmon2","palegreen3","#EEE8AA","cadetblue1","goldenrod1")

p1 <- ggplot(Risk_2021,aes(location,weight = val, fill = rei))+
  geom_bar(color = 'black',width = .7,position = 'dodge',
           size = .2)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = color)+
  theme_classic()+
  coord_flip() + facet_grid(.~rei) + theme_light() +
  geom_text(aes(label=val2, y=val+5.5), 
            position=position_dodge(0.8), vjust=0,hjust = 1,
            size = 2) +
  theme(
        panel.grid = element_blank(),
        panel.spacing.x = unit(0, "cm"),  
        panel.spacing.y = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(family = "Times", color = "black")
        ) +  
  xlab("")+
  ylab("")  
ggsave("global.risk.both.pdf",p1, dpi = 300, width = 14, height = 8, units = "in")



rm(list = ls()) 
library(easyGBDR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
#ls("package:easyGBDR")
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/7.risk') ##设置工作路径

order1 <- read.csv('order.csv',header = F)
order1$V1 <- rev(order1$V1) ### 读取纵坐标不同区域的排列顺序

#?GBDread
data <- GBDread(foldername = '.', #数据所在路径
                 folder = T) %>% 
  filter(measure=="DALYs (Disability-Adjusted Life Years)")
order <- read.csv('order.csv',header = F)
order$V1 <- rev(order$V1) ### 读取纵坐标不同区域的排列顺序

Risk_2021 <- subset(data, year==2021 &
                      metric=='Percent' &
                      sex=='Male' &
                      age=='All ages' )
Risk_2021$val <- round(Risk_2021$val*100,1)
Risk_2021$val2 <- paste0(Risk_2021$val,'%')
Risk_2021$location <- factor(Risk_2021$location, 
                             levels=order$V1, 
                             ordered=TRUE)
Risk_2021_1 <- Risk_2021 %>% filter(location=='Global') %>% arrange(desc(measure),desc(val))

factors  <- Risk_2021_1[1:6,]$rei

Risk_2021$rei <- factor(Risk_2021$rei, 
                        levels= factors, 
                        ordered=TRUE)  ## 按自己想要的顺序排好顺序
Risk_2021 <- na.omit(Risk_2021) 
Risk_2021[which(Risk_2021$measure=='DALYs (Disability-Adjusted Life Years)'),]$measure <-"DALYs"


color <- c("pink2","lightskyblue3","plum3","salmon2","palegreen3","#EEE8AA","cadetblue1","goldenrod1")

p1 <- ggplot(Risk_2021,aes(location,weight = val, fill = rei))+
  geom_bar(color = 'black',width = .7,position = 'dodge',
           size = .2)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = color)+
  theme_classic()+
  coord_flip() + facet_grid(.~rei) + theme_light() +
  geom_text(aes(label=val2, y=val+5.5), 
            position=position_dodge(0.8), vjust=0,hjust = 1,
            size = 2) +
  theme(
        panel.grid = element_blank(),
        panel.spacing.x = unit(0, "cm"),  
        panel.spacing.y = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(family = "Times", color = "black")
  ) +  
  xlab("")+
  ylab("")  
ggsave("global.risk.male.pdf",p1, dpi = 300, width = 14, height = 8, units = "in")


#############################Female####################
rm(list = ls()) 
library(easyGBDR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
#ls("package:easyGBDR")
setwd('C:/Users/yuanfang/Desktop/ISbiomarkers/0.GBD/7.risk') ##设置工作路径

order1 <- read.csv('order.csv',header = F)
order1$V1 <- rev(order1$V1) ### 读取纵坐标不同区域的排列顺序

#?GBDread
data <- data <- GBDread(foldername = '.', #数据所在路径
                        folder = T) %>% 
  filter(measure=="DALYs (Disability-Adjusted Life Years)")
order <- read.csv('order.csv',header = F)
order$V1 <- rev(order$V1) ### 读取纵坐标不同区域的排列顺序

Risk_2021 <- subset(data, year==2021 &
                      metric=='Percent' &
                      sex=='Female' &
                      age=='All ages' )
Risk_2021$val <- round(Risk_2021$val*100,1)
Risk_2021$val2 <- paste0(Risk_2021$val,'%')
Risk_2021$location <- factor(Risk_2021$location, 
                             levels=order$V1, 
                             ordered=TRUE)
Risk_2021_1 <- Risk_2021 %>% filter(location=='Global') %>% arrange(desc(measure),desc(val))

factors  <- Risk_2021_1[1:6,]$rei

Risk_2021$rei <- factor(Risk_2021$rei, 
                        levels= factors, 
                        ordered=TRUE)  ## 按自己想要的顺序排好顺序
Risk_2021 <- na.omit(Risk_2021) 
Risk_2021[which(Risk_2021$measure=='DALYs (Disability-Adjusted Life Years)'),]$measure <-"DALYs"


color <- c("pink2","lightskyblue3","plum3","salmon2","palegreen3","#EEE8AA","cadetblue1","goldenrod1")

p1 <- ggplot(Risk_2021,aes(location,weight = val, fill = rei))+
  geom_bar(color = 'black',width = .7,position = 'dodge',
           size = .2)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = color)+
  theme_classic()+
  coord_flip() + facet_grid(.~rei) + theme_light() +
  geom_text(aes(label=val2, y=val+5.5), 
            position=position_dodge(0.8), vjust=0,hjust = 1,
            size = 2) +
  theme(
        panel.spacing.x = unit(0, "cm"),  
        panel.spacing.y = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(family = "Times", color = "black")
  ) +  
  xlab("")+
  ylab("")  
ggsave("global.risk.female.pdf",p1, dpi = 300, width = 14, height = 8, units = "in")
