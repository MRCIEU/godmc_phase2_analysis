library(ggplot2)
library(tidyverse)
library(ggthemes)
library(meffil)
library(dplyr)
library(gridExtra)
library(stringr)
library(data.table)

load("../results/enrichments/snpcontrolsets_selection.rdata")
w<-which(is.na(f.all$snp_cis))
f.all$Category<-as.character(f.all$snp_cis)
f.all$Category[w]<-"no_mqtl"

f.all$Category<-gsub("TRUE","cisonly SNPs",f.all$Category)
f.all$Category<-gsub("FALSE","transonly SNPs",f.all$Category)
f.all$Category<-gsub("ambivalent","cis+trans SNPs",f.all$Category)

f.all$max_log10pval<-f.all$min_pval
w0<-which(f.all$min_pval==0)
mx<-min(f.all$min_pval[-w0],na.rm=T)
f.all$max_log10pval[w0]<-mx
f.all$max_log10pval<--log10(as.numeric(f.all$max_log10pval))

w<-which(f.all$mqtl_clumped=="TRUE")
f.all2<-f.all[w,]

m<-max(-log10(f.all2$Pvalue))

p1<-ggplot() +
geom_density(data=f.all2, aes(abs(max_abs_Effect),group=Category,color=Category)) +
xlab("Max Abs Effect Size") +
ylab("Density")
ggsave(p1,file="./images/effectsizes.pdf",height=6,width=6)

f.all2$trans_min_pval<-"NA"
f.all2$trans_max_abs_Effect<-"NA"

for (i in 1:23){
cat(i,"\n")
snp_cis<-read_tsv(paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans2.txt.gz"))
w<-which(snp_cis$snp%in%f.all2$SNP)
snp_cis<-snp_cis[w,]
m<-match(f.all2$SNP,snp_cis$snp)
f.all2$trans_min_pval<-snp_cis$trans_min_pval[m]
f.all2$trans_max_abs_Effect<-snp_cis$trans_max_abs_Effect[m]

}

snp_cis<-unique(data.frame(snp=f.all2$SNP,snp_cis=f.all2$Category))
table(snp_cis$snp_cis)

#cisonly SNPs cis+trans SNPs transonly SNPs 
#        157095          66759            794 

##
p1<-ggplot() +
geom_density(data=f.all2, aes(abs(max_log10pval),group=Category,color=Category)) +
xlab("Max -log10 pvals") +
ylab("Density")
ggsave(p1,file="./images/max_log10pvals.pdf",height=6,width=6)

p1<-ggplot() +
geom_histogram(data=f.all2, aes(max_log10pval)) +
facet_wrap(~Category,scale="free_y",ncol=1)+
xlab("Max -log10 pvals")
ggsave(p1,file="./images/max_log10pvals.pdf",height=6,width=6)

####

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped2 <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
table(clumped2$cis)
# FALSE   TRUE 
# 23117 248607 
length(unique(clumped2$cpg))
#190102

data=as.data.table(clumped2)
data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]
data[,Category:=ifelse(all(cis),"cisonly cpgs",ifelse(all(!cis),"trans only cpgs","cis + trans cpgs")),by=c("cpgchr","cpgpos")]

#table(data$Category)

#    cisonly cpgs cis + trans cpgs  trans only cpgs 
#          229924            33158             8642 

cpg_cis<-unique(data.frame(cpg=data$cpg,cpg_cis=data$Category))
table(cpg_cis$cpg_cis)

#    cisonly cpgs cis + trans cpgs  trans only cpgs 
#          170986            11902             7214 

df3<-data[ , (min_Effect = min(Effect)), by = cpg]
df4<-data[ , (max_Effect = max(Effect)), by = cpg]
df<-data.table(rbind(df3,df4))

maxAbsObs <- function(x) x[which.max(abs(x))]
df2<-df[, lapply(.SD, maxAbsObs), by="cpg"]
data<-inner_join(data,df2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"max_abs_Effect"

p1<-ggplot() +
geom_density(data=data, aes(abs(max_abs_Effect),group=Category,color=Category)) +
xlab("Max Abs Effect Size") +
ylab("Density")
ggsave(p1,file="./images/effectsizescpg.pdf",height=6,width=6)

data=as.data.table(data)

minpval<-data[ , (minp = min(pval)), by = cpg]
data<-inner_join(data,minpval)
w<-which(names(data)%in%"V1")
names(data)[w]<-"min_pval"


p1<-ggplot() +
geom_histogram(data=data, aes(-log10(min_pval))) +
facet_wrap(~Category,scale="free_y",ncol=1)+
xlab("Max -log10 pval")
ggsave(p1,file="./images/max_log10pvalscpgs.pdf",height=6,width=6)
