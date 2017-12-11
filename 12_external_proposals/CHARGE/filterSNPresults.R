# Get cpg positions
library(tidyverse)

cis_radius<-as.numeric(1000000)

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/04_conditional_16/cpg_pos.rdata")
#r<-read.table("GlycaemicTraits_SNPs_38cohorts.txt",sep="\t",he=T)

r2.all<-data.frame()
for (CHR in 1:23){
cat(CHR,"\n")
f<-file.info(paste0("chr",CHR))$size

if(f>0){

r<-read.table(paste0("chr",CHR),sep="\t",he=F)
h<-read.table("header.txt",he=T)
names(r)<-names(h)
nrow(r) #23970

r <- r %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
r$snp2 <- r$snp
r <- r %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
r$snppos <- as.numeric(r$snppos)
r <- inner_join(r, cpgpos, by=c("cpg"))
r$cis <- FALSE
r$cis[r$snpchr == r$cpgchr & (abs(r$snppos - r$cpgpos) <= cis_radius)] <- TRUE

s<-read.table("SNPs_GlycaemicTraits.formatted.txt",sep="\t",he=F)
nrow(s) #980

r<-r[r[,1]%in%s[,1],]
length(unique(r$snp)) #804

flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(r$snp%in%flip[,1])
if(length(w)>0){
r<-r[-w,]}

indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
w<-which(r$snp%in%indels[,1])
if(length(w)>0){
r<-r[-w,]}

dim(r) #23970

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509

r<-r[which(r$cpg%in%retaincpg),]
nrow(r)
#[1] 19219

length(unique(r$snp))
#[1] 792

r2<-r[which(r$cis==TRUE),]
dim(r2) #18040
length(unique(r2$snp))
#[1] 792
r2.all<-rbind(r2.all,r2)

}}
write.table(r2.all,"GlycaemicTraits_SNPs_cis_38cohorts17.txt",sep="\t",col.names=T,row.names=F,quote=F)


