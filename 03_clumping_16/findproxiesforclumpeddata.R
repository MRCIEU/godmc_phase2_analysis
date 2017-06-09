library(tidyverse)
library(dplyr)

# Get clumped data
load("../results/16/16_clumped.rdata")
cl<-data.frame(clumped)

cl.out<-data.frame()
for (i in 1:23){
cat(i,"\n")
cl_chr<-cl[which(cl$snpchr%in%paste0("chr",i)),]
a <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.",i,".ld.formatted.gz"))
a<-data.frame(a)

agg<-a %>%  
  group_by(BP_A) %>%  
  summarise (start_bp = min(BP_B),  
             stop_bp = max(BP_B),nproxies = n())

agg<-data.frame(agg)
#a2<-a[which(a$SNP_A%in%cl_chr$snp),]

#test<-a3[a3$BP_A=="16855618",]
#min(test$BP_B)
#[1] 16855677
#max(test$BP_B)
#[1] 16866339

#w<-which(agg$BP_A%in%cl_chr$snppos)
#> length(w)
#[1] 3275

#length(which(is.na(cl_chr$BP_A)))
#[1] 1208

m<-match(cl_chr$snppos,agg$BP_A)
cl_chr<-data.frame(cl_chr,agg[m,])

cl.out<-rbind(cl.out,cl_chr)

}

save(cl.out,file="../results/16/clumpedwithldregion.rdata")


w<-which(is.na(cl.out$BP_A))
length(w)
#[1] 67159
table(miss$snpchr)

# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
# 5903  3666  3907  3435  1730  1804  2013  3222  3827   796  3438  4544  1501 
#chr21 chr22 chr23  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
#  881  1208   269  3093  2795  3462  6651  4278  3297  1439 
