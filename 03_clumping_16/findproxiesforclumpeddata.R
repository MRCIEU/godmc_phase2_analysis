library(tidyverse)
library(dplyr)

# Get clumped data
load("../results/16/16_clumped.rdata")
cl<-data.frame(clumped)

i=22
a <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.",i,".ld.formatted.gz"))
a<-data.frame(a)

a2<-a[which(a$SNP_A%in%cl$snp),]

a2 %>% 
  group_by(as.factor(a2$BP_A)) %>% 
  summarise (min_bp = min(a2$BP_B), 
  	         max_bp = max(a2$BP_B))

#test<-a3[a3$BP_A=="16855618",]
#min(test$BP_B)
#[1] 16855677
#max(test$BP_B)
#[1] 16866339

#library(data.table)
#a3<-a2[1:1000,]
#dt <- data.table(a3)
#dt[,list(min=min(a3$BP_B),max=max(a3$BP_B)),by=a3$BP_A]
#    a3$BP_A      min      max
#1: 16855618 16855677 17065023
#2: 16868784 16855677 17065023
#3: 17057138 16855677 17065023
