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

w<-which(is.na(cl.out$start_bp))
cl.out$start_bp[w]<-cl.out$snppos[w]
cl.out$stop_bp[w]<-cl.out$snppos[w]

w<-which(cl.out$start_bp>cl.out$snppos)
cl.out$start_bp[w]<-cl.out$snppos[w]

w<-which(cl.out$stop_bp<cl.out$snppos)
cl.out$stop_bp[w]<-cl.out$snppos[w]

save(cl.out,file="../results/16/clumpedwithldregion.rdata")


