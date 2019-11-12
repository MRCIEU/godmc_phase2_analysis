arguments<-commandArgs(T)
i<-as.numeric(arguments[1])

library(data.table)
library(dplyr)
###
chr<-paste("chr",i,sep="")
chunks<-c(1:962)

res.chr<-data.frame()
for (chunk in 1:length(chunks)){

load(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",chunk,".rdata",sep=""))
w<-which(res$snpchr==chr)

if(length(w)>0){
res<-res[w,]
res$chunk<-chunk
cat(nrow(res),"\n")
res.chr<-rbind(res.chr,res)}

}

res.chr$id<-paste(res.chr$snp,res.chr$cpg,sep="_")

res.chr<-res.chr[,c("cpg","snp","snpchr","snppos","snptype","id","Pvalue","chunk","cpgchr","cpgpos","cis","Effect","StdErr","TotalSampleSize","Freq1","HetISq")]
#cpg	snp	snpchr	snppos	type	id	pval	chunk	cpgchr	cpgpos	cis	trans	snp_cis	min_pval

#write.table(res.chr,paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".txt",sep=""), sep=" ",col.names=F,row.names=F,quote=F)
#df.out<-read.table(paste("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".txt.gz",sep=""),he=F)
df.out<-res.chr
df.out$cpgchr<-gsub("chrX","chr23",df.out$cpgchr)
names(df.out)<-c("cpg","snp","snpchr","snppos","type","id","pval","chunk","cpgchr","cpgpos","cis","Effect","StdErr","TotalSampleSize","Freq1","HetISq")
#spl<-strsplit(as.character(df.out$id),split="_")
#spl<-do.call("rbind",spl)
#df.out<-data.frame(cpg=spl[,2],snp=spl[,1],df.out)

#retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 #exclusion probes from TwinsUK
#excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
#rm<-which(retaincpg%in%excl[,1])
#14882
#retaincpg<-retaincpg[-rm]
#420509
 
#df.out<-df.out[which(df.out$cpg%in%retaincpg),]

#remove snps
#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(df.out$snp%in%flip[,1])
#if(length(w)>0){
#df.out<-df.out[-w,]}

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(df.out$snp%in%indels[,1]) #129
#if(length(w)>0){
#df.out<-df.out[-w,]}

#library(meffil)
#y<-meffil.get.features("450k")
#m<-match(df.out$cpg,y$name)

#df.out<-data.frame(df.out,y[m,c("chromosome","position")])
#names(df.out)<-c("cpg","snp","snpchr","snppos","type","id","pval","chunk","cpgchr","cpgpos")

#w<-which(df.out$snpchr==df.out$cpgchr&abs(df.out$snppos-df.out$cpgpos)<=1000000)
#df.out$cis2<-"FALSE"
#df.out$cis2[w]<-"TRUE"
#df.out$trans<-"FALSE"
#w<-which(df.out$snpchr==df.out$cpgchr&abs(df.out$snppos-df.out$cpgpos)>1000000)
#df.out$trans[w]<-"TRUE"
#w<-which(df.out$snpchr!=df.out$cpgchr)
#df.out$trans[w]<-"TRUE"

#
w1<-which(df.out$trans==TRUE & df.out$pval > 1e-14)
w2<-which(df.out$cis==TRUE & df.out$pval > 1e-8)
dim(df.out)
rm<-unique(c(w1,w2))
df.out<-df.out[-rm,]
dim(df.out)

#test
#w<-which(df.out$snp=="chr22:49697538:SNP")

#              cpg                snp snpchr   snppos type
#6      cg00021237 chr22:49697538:SNP  chr22 49697538  SNP
#337076 cg06629999 chr22:49697538:SNP  chr22 49697538  SNP
#349035 cg06847006 chr22:49697538:SNP  chr22 49697538  SNP
#802869 cg19521610 chr22:49697538:SNP  chr22 49697538  SNP
#                                  id       pval chunk cpgchr   cpgpos   cis
#6      chr22:49697538:SNP_cg00021237 5.370e-236     1   chr8 67525849 FALSE
#337076 chr22:49697538:SNP_cg06629999  1.071e-21   255  chr22 49720192  TRUE
#349035 chr22:49697538:SNP_cg06847006 1.920e-134   263  chr22 49697565  TRUE
#802869 chr22:49697538:SNP_cg19521610  1.059e-34   701  chr22 49717457  TRUE
 
data<-data.table(df.out)
data$cis<-as.logical(data$cis)
df1<-data[,snp_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("snp")]

#data[w,]
#          cpg                snp snpchr   snppos type
#1: cg00021237 chr22:49697538:SNP  chr22 49697538  SNP
#2: cg06629999 chr22:49697538:SNP  chr22 49697538  SNP
#3: cg06847006 chr22:49697538:SNP  chr22 49697538  SNP
#4: cg19521610 chr22:49697538:SNP  chr22 49697538  SNP
#                              id       pval chunk cpgchr   cpgpos   cis trans
#1: chr22:49697538:SNP_cg00021237 5.370e-236     1   chr8 67525849 FALSE  TRUE
#2: chr22:49697538:SNP_cg06629999  1.071e-21   255  chr22 49720192  TRUE FALSE
#3: chr22:49697538:SNP_cg06847006 1.920e-134   263  chr22 49697565  TRUE FALSE
#4: chr22:49697538:SNP_cg19521610  1.059e-34   701  chr22 49717457  TRUE FALSE
#      snp_cis
#1: ambivalent
#2: ambivalent
#3: ambivalent
#4: ambivalent

data[which(data$snp%in%c("chr22:29654004:SNP")),]

df2<-data[ , (min_pval = min(pval)), by = snp]
data<-inner_join(df1,df2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"min_pval"

#data<-data.table(data)
#df3<-data[ , (min_Effect = min(Effect)), by = snp]
#df4<-data[ , (max_Effect = max(Effect)), by = snp]
#df<-data.table(rbind(df3,df4))

data<-data %>% 
  group_by(snp) %>% 
  slice(which.max(abs(Effect)))

w<-which(names(data)%in%"Effect")
names(data)[w]<-"max_abs_Effect"

w<-which(names(data)%in%"StdErr")
names(data)[w]<-"max_abs_SE"

#maxAbsObs <- function(x) x[which.max(abs(x))]
#df2<-df[, lapply(.SD, maxAbsObs), by="snp"]
#data<-inner_join(data,df2)
#w<-which(names(data)%in%"V1")
#names(data)[w]<-"max_abs_Effect"

write.table(data,paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans3.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
#
###
#data<-data.table(data)
#amb<-subset(data,snp_cis!="TRUE" & cis!="TRUE")
#df2<-amb[ , (min_transpval = min(pval)), by = snp]
#data<-full_join(data,df2)
#w<-which(names(data)%in%"V1")
#names(data)[w]<-"trans_min_pval"

#df3<-amb[ , (min_Effect = min(Effect)), by = snp]
#df4<-amb[ , (max_Effect = max(Effect)), by = snp]
#df<-data.table(rbind(df3,df4))

#maxAbsObs <- function(x) x[which.max(abs(x))]
#df2<-df[, lapply(.SD, maxAbsObs), by="snp"]
#data<-full_join(data,df2)
#w<-which(names(data)%in%"V1")
#names(data)[w]<-"trans_max_abs_Effect"

#chr22:50546868:SNP

#write.table(data,paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans2.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
#cpg	snp	snpchr	snppos	type	id	pval	chunk	cpgchr	cpgpos	cis	Effect	snp_cis	min_pval	max_abs_Effect	trans_min_pval trans_max_abs_Effect
