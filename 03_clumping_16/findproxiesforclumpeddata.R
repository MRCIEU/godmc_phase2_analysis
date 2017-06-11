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

agg1<-a %>%  
  group_by(BP_A) %>%  
  summarise (bpa = first(BP_A),start_bp = min(BP_B),  
             stop_bp = max(BP_B),nproxies = n(),chr=first(CHR_A))

agg2<-a %>%  
  group_by(BP_B) %>%  
  summarise (bpb = first(BP_B),start_bp = min(BP_A),  
             stop_bp = max(BP_A),nproxies = n(),chr=first(CHR_B))

b <- full_join(agg1, agg2, c("BP_A" = "BP_B"))

#b$final_min <- b$start_bp.x
#b$final_max <- b$stop_bp.y
#b$final_min[is.na(b$final_min)] <- b$bpa[is.na(b$final_min)]
#b$final_max[is.na(b$final_max)] <- b$bpb[is.na(b$final_max)]
b<-data.frame(b)

b2<-data.frame(b$nproxies.x,b$nproxies.y)
b$nproxies<-apply(b2,1,function(x) y<-sum(x,na.rm=T))

b2<-data.frame(b$start_bp.x,b$start_bp.y)
b$final_min<-apply(b2,1,function(x) y<-min(x,na.rm=T))

b2<-data.frame(b$stop_bp.x,b$stop_bp.y)
b$final_max<-apply(b2,1,function(x) y<-max(x,na.rm=T))


b$CHR<-b$chr.x
b$CHR[is.na(b$CHR)]<-b$chr.y[is.na(b$CHR)]
b$CHR<-paste0("chr",b$CHR)

table(is.na(b$final_min))
table(is.na(b$final_max))
table(b$final_min <= b$final_max)
table(b$final_min <= b$final_max)


b <- select(b, CHR=CHR,SNP=BP_A, min=final_min, max=final_max,nproxies)
summary(b$max - b$min)


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


