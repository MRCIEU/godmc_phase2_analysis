library(tidyverse)
library(ggthemes)
library(meffil)
library(dplyr)
library(gridExtra)
library(stringr)

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
max(clumped[which(clumped$cis==TRUE),"pval"])
#1e-4
max(clumped[which(clumped$cis==FALSE),"pval"])
#5e-8

flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(clumped$snp%in%flip[,1])
clumped<-clumped[-w,]

w1<-which(clumped$snptype=="INDEL")
w2<-which(clumped$snptype=="SNP")
mean(clumped$HetISq[w1])
#[1] 45.66107
mean(clumped$HetISq[w2])
#[1] [1] 45.39092


indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
w<-which(clumped$snp%in%indels[,1]) #129
indels<-data.frame(clumped[w,])
mean(indels$HetISq)
#[1] 58.1876

clumped<-clumped[-w,]

w1<-which(clumped$snptype=="INDEL")
w2<-which(clumped$snptype=="SNP")
mean(clumped$HetISq[w1])
#[1] 45.58455
mean(clumped$HetISq[w2])
#[1] 45.39092


retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
 
clumped<-clumped[which(clumped$cpg%in%retaincpg),]
nrow(clumped)

r<-read.table("CpGs_GlycaemicTriats.txt",sep="\t",he=F)
w<-which(clumped$cpg%in%r$V1)

subset<-clumped[w,]
subset<-subset[subset$cis==TRUE,]
length(unique(subset$cpg)) #1035
dim(subset) #1447

write.table(subset,"GlycaemicTraits_CpGs_clumped_cis_38cohorts.txt",sep="\t",col.names=T,row.names=F,quote=F)

s<-read.table("SNPs_GlycaemicTraits.txt",sep="\t",he=F)
s[,1]<-gsub("{SNP}","SNP",s[,1],fixed=T)
s<-paste0("chr",s[,1])
write.table(s,"SNPs_GlycaemicTraits.formatted.txt",sep="\t",quote=F,row.names=F,col.names=F)


