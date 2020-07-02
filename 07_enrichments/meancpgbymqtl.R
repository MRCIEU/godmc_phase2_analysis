library(GenomicRanges)
library(data.table)
library(tidyr)
library(meffil)
library(SDMTools)
library(dplyr)
library(gridExtra)
library(grid)

arguments<-commandArgs(T)
i<-as.numeric(arguments[1])
chr<-paste0("chr",i)


##load cell type conversion and colors

#
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
max(clumped[which(clumped$cis==TRUE),"pval"])
#1e-4
max(clumped[which(clumped$cis==FALSE),"pval"])
#5e-8


#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(clumped$snp%in%flip[,1])
#clumped<-clumped[-w,]

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(clumped$snp%in%indels[,1]) #129
#clumped<-clumped[-w,]

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
 
#clumped<-clumped[which(clumped$cpg%in%retaincpg),]
#nrow(clumped)

clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
data=as.data.table(clumped)


data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]
#data[,cpg_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]

####run with internal background
#grs_mqtl_effect=create_grs(data=data,selector=c("cpg_change=='mqtl_effect>0'","cpg_change=='mqtl_effect<0'","cpg_change=='ambivalent'","snp_change=='mqtl_effect>0'","snp_change=='mqtl_effect<0'","snp_change=='ambivalent'"))
#produceLOLA_plots_intern(grs=grs_mqtl_effect,plot_pref="mqtl_effect",height=35,width=18)


data[,cpg_cis:=ifelse(all(cis),"cis only",ifelse(all(!cis),"trans only","cis+trans")),by=c("cpg")]
r<-unique(data.frame(cpg=data$cpg,cpgchr=data$cpgchr,min=as.numeric(data$cpgpos),max=as.numeric(data$cpgpos),cpg_cis=data$cpg_cis))
gr_cpg_a = with(r,GRanges(seqnames=cpgchr,ranges=IRanges(min,max),strand=Rle("+")))

r$cpgchr<-gsub("chrX","chr23",r$cpgchr)
r$cpgchr<-gsub("chrY","chr24",r$cpgchr)

#
gr_cpg = with(r,GRanges(seqnames=cpgchr,ranges=IRanges(min,max),strand=Rle("+")))

#
cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda.txt")
ss$V1<-gsub("_16","",ss$V1)
ss$V1<-gsub("00_ARIES","ARIES",ss$V1)
ss<-ss[which(ss$V1%in%c("MARS_omni","Factor_V_Leiden_Family_Study")==F),]
cohort<-ss$V1

y<-meffil.get.features("450k")
cpgs<-unique(y$name)
sd.out<-data.frame(cpgs)
mean.out<-data.frame(cpgs)

for (i in 1:length(cohort)){
cat(cohort[i],"\n")
load(paste0(cohort_dir,cohort[i],"_01/results/01/methylation_summary.RData"))
m<-match(cpgs,row.names(meth_summary))
sd<-meth_summary[m,"sd"]
sd.out<-data.frame(sd.out,sd)
mean.df<-meth_summary[m,"mean"]
mean.out<-data.frame(mean.out,mean.df)

}
#row.names(sd.out)<-cpgs
row.names(mean.out)<-cpgs

names(sd.out)[-1]<-cohort
sd.out<-sd.out[,-1]

names(mean.out)[-1]<-cohort
mean.out<-mean.out[,-1]


#mean.out<-as.matrix(mean.out)
#wmean<-rowWeightedMeans(mean.out,w=ss$V2)

#sdmean<-data.frame(cpgs,sd=rowMeans(sd.out, na.rm=TRUE))
#meanmean<-data.frame(cpgs,mean=rowMeans(mean.out, na.rm=TRUE))

#pool_mean_sd <- function(n, m, s)
#{
#    gn <- sum(n, na.rm=TRUE)
#    gm <- sum(m * n, na.rm=TRUE) / gn
#    gs <- sqrt(sum(n * (s^2 + (m - gm)^2), na.rm=TRUE) / gn)
#    return(data.frame(gn=gn, gm=gm, gs=gs))
#}

#p.sd<-pool_mean_sd(n=ss$V2,m=mean.out[1,],s=sd.out[1,])
#a<-apply(mean.out,1,function(x) pool_mean_sd(n=ss$V2,m=x,s=sd.out))

meanmean<-apply(mean.out,1,function(x) wt.mean(as.numeric(x),wt=as.numeric(ss$V2)))
sdmean<-apply(sd.out,1,function(x) wt.sd(as.numeric(x),wt=as.numeric(ss$V2)))

m<-match(data$cpg,names(meanmean))
df2<-unique(data.frame(cpg=data$cpg,meancpg=meanmean[m],sdcpg=sdmean[m],cpg_cis=data$cpg_cis))

m<-match(names(meanmean),data$cpg)
df.all<-unique(data.frame(cpg=names(meanmean),meancpg=meanmean,sdcpg=sdmean,cpg_cis=data$cpg_cis[m]))
m<-match(retaincpg,df.all$cpg)
df.all<-df.all[m,]
save(df.all,file="mean_allcpgs.Robj")

df.all %>%
  group_by(cpg_cis) %>%
  dplyr::summarize(SD = mean(sdcpg, na.rm=TRUE),mean = mean(meancpg, na.rm=TRUE),median=median(meancpg,na.rm=TRUE))

# A tibble: 4 x 3
#cpg_cis         SD      mean
#<fctr>      <dbl>     <dbl>
#  1    no mqtl 0.01254571 0.5339780
#2   cis only 0.01344013 0.4715994
#3  cis+trans 0.01081446 0.4411792
#4 trans only 0.01036164 0.3452771
df3 %>%
  group_by(cpg_cis) %>%
  dplyr::summarize(SD = median(sdcpg, na.rm=TRUE),mean = mean(meancpg, na.rm=TRUE),median=median(meancpg,na.rm=TRUE))


data$MAF<-data$Freq1
w<-which(data$MAF>0.5)
data$MAF[w]<-1-data$MAF[w]


###
data2<-data[ , (min_pval = min(pval)), by = cpg]
data<-inner_join(data,data2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"min_pval"

data<-data.table(data)
data2<-data[ , (min_maf = min(MAF)), by = cpg]
data<-inner_join(data,data2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"min_maf"

data<-data.table(data)
data2<-data[ , (max_maf = min(MAF)), by = cpg]
data<-inner_join(data,data2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"max_maf"

data<-data.table(data)
data2<-data[ , (max_i2 = max(HetISq)), by = cpg]
data<-inner_join(data,data2)
w<-which(names(data)%in%"V1")
names(data)[w]<-"max_i2"

data<-data.table(data)
df3<-data[ , (min_Effect = min(Effect)), by = cpg]
df4<-data[ , (max_Effect = max(Effect)), by = cpg]
df3<-data.table(rbind(df3,df4))

maxAbsObs <- function(x) x[which.max(abs(x))]
df3<-df3[, lapply(.SD, maxAbsObs), by="cpg"]
data<-inner_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"max_abs_Effect"
###
data<-data.table(data)

amb<-subset(data,cpg_cis!="cis only" & cis!="TRUE")
df3<-amb[ , (min_transpval = min(pval)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"trans_min_pval"

df3<-amb[ , (min_transmaf = max(MAF)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"trans_max_maf"

df3<-amb[ , (min_cismaf = min(MAF)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"trans_min_maf"

df3<-amb[ , (max_transi2 = max(HetISq)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"trans_max_i2"

df3<-amb[ , (min_Effect = min(Effect)), by = cpg]
df4<-amb[ , (max_Effect = max(Effect)), by = cpg]
df3<-data.table(rbind(df3,df4))

maxAbsObs <- function(x) x[which.max(abs(x))]
df3<-df3[, lapply(.SD, maxAbsObs), by="cpg"]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"trans_max_abs_Effect"


##
data<-data.table(data)

amb<-subset(data,cpg_cis!="trans only" & cis!="FALSE")
df3<-amb[ , (min_cispval = min(pval)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"cis_min_pval"

df3<-amb[ , (min_cismaf = max(MAF)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"cis_max_maf"

df3<-amb[ , (min_cismaf = min(MAF)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"cis_min_maf"

df3<-amb[ , (max_cisi2 = max(HetISq)), by = cpg]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"cis_max_i2"

df3<-amb[ , (min_Effect = min(Effect)), by = cpg]
df4<-amb[ , (max_Effect = max(Effect)), by = cpg]
df3<-data.table(rbind(df3,df4))

maxAbsObs <- function(x) x[which.max(abs(x))]
df3<-df3[, lapply(.SD, maxAbsObs), by="cpg"]
data<-full_join(data,df3)
w<-which(names(data)%in%"V1")
names(data)[w]<-"cis_max_abs_Effect"
##

tmp<-unique(data[,c("cpg","min_pval","max_abs_Effect","trans_min_pval","trans_max_abs_Effect","cis_min_pval","cis_max_abs_Effect","min_maf","cis_min_maf","trans_min_maf","max_i2","cis_max_i2","trans_max_i2")])
m<-match(df2$cpg,tmp$cpg)
df2<-data.frame(df2,tmp[m,-1])

m<-match(data$cpg,df2$cpg)
data<-data.frame(data,df2[m,c("meancpg","sdcpg")])

df<-data.frame(df2[,1:3],cpg_cis="All",df2[,5:ncol(df2)])
df2<-rbind(df,df2)



amb_cis<-df2[df2$cpg_cis=="cis+trans",]
amb_trans<-df2[df2$cpg_cis=="cis+trans",]
amb_cis$cpg_cis<-"cis+trans_cis"
amb_trans$cpg_cis<-"cis+trans_trans"

amb_cis$max_abs_Effect<-amb_cis$cis_max_abs_Effect
amb_trans$max_abs_Effect<-amb_trans$trans_max_abs_Effect
amb_cis$min_maf<-amb_cis$cis_min_maf
amb_trans$min_maf<-amb_trans$trans_min_maf
amb_cis$max_i2<-amb_cis$cis_max_i2
amb_trans$max_i2<-amb_trans$trans_max_i2


df3<-rbind(df2,amb_cis,amb_trans)
o<-order(df3$cpg_cis)
df3<-df3[o,]
w<-which(df3$cpg_cis%in%c("cis+trans"))
df3<-df3[-w,]
#p1 <- ggplot(df2, aes(x=as.factor(cpg_cis), y=meancpg)) +
#geom_boxplot() +
#theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
#labs(x="Category", y="mean methylation")
#ggsave(plot=p1, file="./images/meancpg_ciscategory.pdf", width=7, height=7)

m<-match(df2$cpg,y$name)
df2<-data.frame(df2,y[m,])

n<-data.frame(table(df2$cpg_cis))
df2$cpg_cis_n<-df2$cpg_cis
df2$cpg_cis_n<-gsub("All",paste0("All (n=",n$Freq[1],")"),df2$cpg_cis_n)
df2$cpg_cis_n<-gsub("cis only",paste0("cis only (n=",n$Freq[2],")"),df2$cpg_cis_n)
df2$cpg_cis_n<-gsub("cis\\+trans",paste0("cis+trans (n=",n$Freq[3],")"),df2$cpg_cis_n)
df2$cpg_cis_n<-gsub("trans only",paste0("trans only (n=",n$Freq[4],")"),df2$cpg_cis_n)
save(df2,file="meancpg.Robj")

m<-match(df.all$cpg,y$name)
df.all<-data.frame(df.all,y[m,])

w<-which(is.na(df.all$cpg_cis))
df.all$cpg_cis<-as.character(df.all$cpg_cis)
df.all$cpg_cis[w]<-"No mqtl"
df.all$cpg_cis<-as.factor(df.all$cpg_cis)
n<-data.frame(table(df.all$cpg_cis))
df.all$cpg_cis_n<-df.all$cpg_cis
df.all$cpg_cis_n<-gsub("No mqtl",paste0("no mQTL (n=",n$Freq[3],")"),df.all$cpg_cis_n)
df.all$cpg_cis_n<-gsub("cis only",paste0("cis only (n=",n$Freq[1],")"),df.all$cpg_cis_n)
df.all$cpg_cis_n<-gsub("cis\\+trans",paste0("cis+trans (n=",n$Freq[2],")"),df.all$cpg_cis_n)
df.all$cpg_cis_n<-gsub("trans only",paste0("trans only (n=",n$Freq[4],")"),df.all$cpg_cis_n)

p1 <- ggplot(df2, aes(x=meancpg)) +
geom_histogram() +
facet_wrap(~cpg_cis_n,scales="free_y") +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),strip.text = element_text(size = 8)) +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean.pdf", width=7, height=7)

p1<-ggplot(df2,aes(x=meancpg,fill=relation.to.island)) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
facet_wrap(~cpg_cis_n,scales="free_y") +
ylab("Proportion CpGs") +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean_cpgisland.pdf", width=7, height=7)


df.all$facet = factor(df.all$cpg_cis_n, levels = c("no mQTL (n=230407)", "cis only (n=170986)", "cis+trans (n=11902)", "trans only (n=7214)"))

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(org.Hs.eg.db)
library(limma)
library(missMethyl)
library(scales)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#m<-match(df.all$cpg,ann$Name)
#df.all<-data.frame(df.all,ann[,"Regulatory_Feature_Group"])
df.all<-data.frame(merge(df.all,ann[,c(4,32)],by.x="cpg",by.y="Name",all.x=T))

load("cpgbysegstate.Robj")
m<-match(st$cpg,df.all$cpg)
df.all2<-data.frame(df.all[m,],st)
w<-which(is.na(df.all2$facet))
df.all2<-df.all2[-w,]

p1<-ggplot(df.all2,aes(x=meancpg,fill=state)) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
facet_wrap(~facet,scales="free_y",nr=2) +
ylab("Proportion CpGs") +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean_nomqtl_state.pdf", width=7, height=7)

p1<-ggplot(df.all2, aes(fill=state, y=1, x=facet)) + 
geom_bar( stat="identity", position="fill") +
labs(x="Category", y="proportion CpGs")
ggsave(plot=p1, file="./images/statecpg_ciscategory.pdf", width=7, height=7)

p1 <- ggplot(df7, aes(x=as.factor(cpg_cis), y=1,fill=gene.annotation)) +
geom_bar(stat="identity") +
scale_y_continuous(labels = percent_format())
ggsave(plot=p1, file="./images/geneannotationcpg_ciscategory.pdf", width=7, height=7)

p1<-ggplot(df7, aes(fill=gene.annotation, y=1, x=facet)) + 
geom_bar( stat="identity", position="fill") +
labs(x="Category", y="proportion CpGs")
ggsave(plot=p1, file="./images/geneannotationcpg_ciscategory.pdf", width=7, height=7)

p1<-ggplot(df.all,aes(x=meancpg,fill=relation.to.island)) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
#facet_wrap(~facet,scales="free_y",nr=1) +
facet_wrap(~facet,nr=1) +
labs(y="Proportion DNAm sites",x="Weighted mean by mQTL site",fill = "DNAm site annotation") +
theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14),axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18), strip.text = element_text(size = 20)) +
theme(legend.text=element_text(size=20),legend.title=element_text(size=20)) +
scale_fill_brewer(type="qual")
ggsave(plot=p1, file="./images/cpgmean_nomqtl_cpgisland.pdf", width=18, height=6)

p1<-ggplot(df.all,aes(x=meancpg,fill=type)) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
facet_wrap(~facet,scales="free_y",nr=2) +
ylab("Proportion CpGs") +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean_nomqtl_typei_ii.pdf", width=7, height=7)

p1<-ggplot(df.all,aes(x=meancpg,fill=Regulatory_Feature_Group)) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
facet_wrap(~facet,scales="free_y",nr=2) +
ylab("Proportion CpGs") +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean_nomqtl_regulatoryfeature.pdf", width=7, height=7)


g1<-grep("5'UTR",df.all$gene.region)
g2<-grep("1stExon",df.all$gene.region)
g3<-grep("TSS200",df.all$gene.region)
g4<-grep("TSS1500",df.all$gene.region)
g5<-grep("3'UTR",df.all$gene.region)
g6<-grep("Body",df.all$gene.region)

g7<-unique(c(g1,g2,g3,g4,g5,g6))
g7<-data.frame(df.all[-g7,],gene.annotation="Intergenic")

g1<-data.frame(df.all[grep("5'UTR",df.all$gene.region),],gene.annotation="5'UTR")
g2<-data.frame(df.all[grep("1stExon",df.all$gene.region),],gene.annotation="1stExon")
g3<-data.frame(df.all[grep("TSS200",df.all$gene.region),],gene.annotation="TSS200")
g4<-data.frame(df.all[grep("TSS1500",df.all$gene.region),],gene.annotation="TSS1500")
g5<-data.frame(df.all[grep("3'UTR",df.all$gene.region),],gene.annotation="3'UTR")
g6<-data.frame(df.all[grep("Body",df.all$gene.region),],gene.annotation="Body")

df7<-rbind(g1,g2,g3,g4,g5,g6,g7)

p1<-ggplot(df7, aes(x = meancpg,fill=gene.annotation)) + 
#geom_histogram(aes(y=..count../sum(..count..))) + 
geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
facet_wrap(~facet,scales="free_y") +
ylab("Proportion CpGs") +
xlab("weighted mean by cpg")
ggsave(plot=p1, file="./images/cpgmean_nomqtl_generegion.pdf", width=7, height=7)

p1 <- ggplot(df7, aes(x=as.factor(cpg_cis), y=1,fill=gene.annotation)) +
geom_bar(stat="identity") +
scale_y_continuous(labels = percent_format())
ggsave(plot=p1, file="./images/geneannotationcpg_ciscategory.pdf", width=7, height=7)

p1<-ggplot(df7, aes(fill=gene.annotation, y=1, x=facet)) + 
geom_bar( stat="identity", position="fill") +
labs(x="Category", y="proportion CpGs")
ggsave(plot=p1, file="./images/geneannotationcpg_ciscategory.pdf", width=7, height=7)

table(df2$cpg_cis)
#All   cis only  cis+trans trans only
#190102     170986      11902       7214

table(df3$cpg_cis)
#            All        cis only       cis+trans      trans only   cis+trans_cis 
#         190102          170986           11902            7214           11902 
#cis+trans_trans 
#          11902 

table(df.all$relation.to.island,df.all$cpg_cis)
#cis only cis+trans No mqtl trans only
#                 0         0      65          0
#  Island     37287      3138   91635       3686
#  N_Shelf     8981       508   10659        267
#  N_Shore    26500      1887   26550        716
#  OpenSea    69121      4444   71339       1762
#  S_Shelf     8105       382    9566        203
#  S_Shore    20992      1543   20593        580



p1 <- ggplot(df3, aes(x=as.factor(cpg_cis), y=min_maf)) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="min MAF")
ggsave(plot=p1, file="./images/minmaf_cpg_ciscategory.pdf", width=7, height=7)

p1 <- ggplot(df3, aes(x=as.factor(cpg_cis), y=max_i2)) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="max I2")
ggsave(plot=p1, file="./images/maxi2_cpg_ciscategory.pdf", width=7, height=7)

df3$cpg_cis<-gsub("All","all",df3$cpg_cis)
p1 <- ggplot(df3, aes(x=as.factor(cpg_cis), y=abs(max_abs_Effect),fill=as.factor(cpg_cis))) +
#geom_boxplot() +
geom_violin() +
geom_boxplot(width=0.1,fill="white") +
theme(axis.title.y=element_text(size=20),axis.title.x=element_text(size=20),axis.text.x=element_text(angle = 90, hjust = 1,size=16),axis.text.y=element_text(size=16)) +
labs(x="mQTL category", y="Maximum absolute mQTL effect size") +
theme(legend.position="none")
ggsave(plot=p1, file="./images/effectsizescpg_ciscategory.pdf", width=7, height=7)

summary(lm(abs(df3[w,"max_abs_Effect"])~df3[w,"cpg_cis"]))

#Coefficients:
#                              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                  0.2850567  0.0005418 526.145  < 2e-16 ***
#df3[w, "cpg_cis"]trans only -0.0205893  0.0026927  -7.646 2.08e-14 ***

#table(df3[w, "cpg_cis"])

#            All        cis only       cis+trans      trans only   cis+trans_cis 
#              0          170986               0            7214               0 
#cis+trans_trans 
#              0 

df.all$cpg_cis<-gsub("No mqtl","no mQTL",df.all$cpg_cis)
df.all$cpg_cis<-gsub("cis+trans      ","cis+trans",df.all$cpg_cis)
df.all$cpg_cis = factor(df.all$cpg_cis, levels = c("no mQTL", "cis only", "cis+trans", "trans only"))

p2 <- ggplot(df.all, aes(x=as.factor(cpg_cis), y=sdcpg,fill=as.factor(cpg_cis))) +
geom_violin() +
geom_boxplot(width=0.1,fill="white") +
theme(axis.title.y=element_text(size=20),axis.title.x=element_text(size=20),axis.text.x=element_text(angle = 90, hjust = 0,size=16),axis.text.y=element_text(size=16)) +
labs(x="mQTL category", y="Weighted SD methylation") +
theme(legend.position="none")
ggsave(plot=p2, file="./images/sdcpg_ciscategory.pdf", width=7, height=7)

pdf("FigSX_effectsize_SD.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
dev.off()

p1 <- ggplot(df2, aes(x=as.factor(cpg_cis), y=abs(trans_max_abs_Effect))) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="max abs trans Effect")
ggsave(plot=p1, file="./images/transeffectsizescpg_ciscategory.pdf", width=7, height=7)

p1 <- ggplot(df2, aes(x=as.factor(cpg_cis), y=abs(cis_max_abs_Effect))) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="max abs cis Effect")
ggsave(plot=p1, file="./images/ciseffectsizescpg_ciscategory.pdf", width=7, height=7)

p1 <- ggplot(df2, aes(x=sdcpg)) +
geom_histogram() +
facet_wrap(~cpg_cis,scales="free_y") +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs( x="mean cpg")
ggsave(plot=p1, file="./images/sdmean.pdf", width=7, height=7)

p1 <- ggplot(df2, aes(x=abs(cis_max_abs_Effect))) +
geom_histogram() +
facet_wrap(~cpg_cis,scales="free_y") +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs( x="mean cpg")
ggsave(plot=p1, file="./images/ciseffectsizescpg_ciscategory.pdf", width=7, height=7)

p1 <- ggplot(df2, aes(x=abs(trans_max_abs_Effect))) +
geom_histogram() +
facet_wrap(~cpg_cis,scales="free_y") +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs( x="mean cpg")
ggsave(plot=p1, file="./images/transeffectsizescpg_ciscategory.pdf", width=7, height=7)


df5<-unique(data.frame(cpg=data$cpg,min_maf=data$min_maf,cpg_cis=data$cpg_cis,cis=data$cis))
df6<-df5[which(df5$cpg$cpg_cis=="cis+trans"),]



df2%>%group_by(cpg_cis)%>%summarize(mean=mean(meancpg))
#A tibble: 4 x 2
#     cpg_cis      mean
#      <fctr>     <dbl>
#1        All 0.5229117
#2   cis only 0.5339780
#3  cis+trans 0.4715994
#4 trans only 0.3452771

df2%>%group_by(cpg_cis)%>%summarize(sd=mean(sdcpg))
# A tibble: 4 x 2
#     cpg_cis         sd
#      <fctr>      <dbl>
#1        All 0.01251882
#2   cis only 0.01254571
#3  cis+trans 0.01344013
#4 trans only 0.01036164

df2%>%group_by(cpg_cis)%>%summarize(sd=median(sdcpg))

df2%>%group_by(cpg_cis)%>%summarize(es=median(abs(trans_max_abs_Effect)))
# A tibble: 4 x 2
#     cpg_cis        es
#      <fctr>     <dbl>
#1        All        NA
#2   cis only        NA
#3  cis+trans 0.1701506
#4 trans only 0.2010279
df2%>%group_by(cpg_cis)%>%summarize(es=median(abs(cis_max_abs_Effect)))
# A tibble: 4 x 2
#     cpg_cis        es
#      <fctr>     <dbl>
#1        All        NA
#2   cis only 0.2115842
#3  cis+trans 0.2708918
#4 trans only        NA
df2%>%group_by(cpg_cis)%>%summarize(es=median(abs(max_abs_Effect)))
# A tibble: 4 x 2
#     cpg_cis        es
#      <fctr>     <dbl>
#1        All 0.2178641
#2   cis only 0.2115842
#3  cis+trans 0.3181236
#4 trans only 0.2010279

df2%>%group_by(cpg_cis)%>%summarize(es=min(abs(max_abs_Effect)))
# A tibble: 4 x 2
 #    cpg_cis        es
 #     <fctr>     <dbl>
#1        All 0.0423828
#2   cis only 0.0423828
#3  cis+trans 0.0741861
#4 trans only 0.0456189

df2%>%group_by(cpg_cis)%>%summarize(es=max(abs(max_abs_Effect)))
# A tibble: 4 x 2
#     cpg_cis       es
#      <fctr>    <dbl>
#1        All 2.157811
#2   cis only 2.157811
#3  cis+trans 1.799202
#4 trans only 1.950732

table(df2$cpg_cis)
#       All   cis only  cis+trans trans only 
#    190102     170986      11902       7214 

save(data,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped_cpgciscategories.rdata")

summary(lm(abs(df2$max_abs_Effect)~df2$cpg_cis))
#Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            0.2900509  0.0005166 561.466  < 2e-16 ***
#df2$cpg_ciscis only   -0.0049942  0.0007507  -6.653 2.88e-11 ***
#df2$cpg_ciscis+trans   0.0872540  0.0021282  40.998  < 2e-16 ***
#df2$cpg_cistrans only -0.0255835  0.0027017  -9.469  < 2e-16 ***


cis<-df2[which(df2$cpg_cis%in%c("cis only","cis+trans")),]
summary(lm(abs(cis$max_abs_Effect)~cis$cpg_cis))

#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.2850567  0.0005465  521.58   <2e-16 ***
#cis$cpg_ciscis+trans 0.0922482  0.0021423   43.06   <2e-16 ***

w<-which(df2$cpg_cis%in%c("cis only","trans only"))
summary(lm(abs(df2$max_abs_Effect[w])~df2$cpg_cis[w]))
#Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
#(Intercept)               0.2850567  0.0005418 526.145  < 2e-16 ***
#df2$cpg_cis[w]trans only -0.0205893  0.0026927  -7.646 2.08e-14 ***

cis<-df3[which(df3$cpg_cis%in%c("cis only","cis+trans_cis")),]
summary(lm(abs(cis$max_abs_Effect)~cis$cpg_cis))
#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
#(Intercept)              0.2850567  0.0005475  520.63   <2e-16 ***
#cis$cpg_ciscis+trans_cis 0.0520097  0.0021463   24.23   <2e-16 ***

trans<-df3[which(df3$cpg_cis%in%c("trans only","cis+trans_trans")),]
summary(lm(abs(trans$max_abs_Effect)~trans$cpg_cis))

#Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             0.209010   0.001443  144.85   <2e-16 ***
#trans$cpg_cistrans only 0.055457   0.002349   23.61   <2e-16 ***



load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/05_cis-trans-networks/results/graph.rdata")
w1<-which(df2$cpg%in%unique(dat$creg))
w2<-which(df2$cpg%in%unique(dat$tcpg))
df2$creg<-"FALSE"
df2$tcpg<-"FALSE"
df2$creg[w1]<-"TRUE"
df2$tcpg[w2]<-"TRUE"

p1 <- ggplot(df2, aes(x=meancpg)) +
geom_histogram() +
facet_wrap(~creg,scales="free_y") +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs( x="mean cpg")
ggsave(plot=p1, file="./images/cpgmeancreg.pdf", width=7, height=4)

p1 <- ggplot(df2, aes(x=meancpg)) +
geom_histogram() +
facet_wrap(~tcpg,scales="free_y") +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs( x="mean cpg")
ggsave(plot=p1, file="./images/cpgmeantcpg.pdf", width=7, height=4)

#collapse overlaps
path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs"
l<-list.files(path)
l2<-gsub("tfbs_","",l)

ov<-data.frame(r$cpg)
n.out<-NULL
for (i in 1:length(l)){
#for (i in 1:10){
cat(i,"\n")
bed<-read.table(paste0(path,"/",l[i]))
bed$V1<-gsub("chrX","chr23",bed$V1)
bed$V1<-gsub("chrY","chr24",bed$V1)
#n<-paste(bed$V4[1],bed$V5[1],sep="_")
#n.out<-append(n.out,n)
bed<-bed[which(bed$V1==chr),]

#     V1    V2    V3   V4    V5 V6
#1 chr21 11529 11731 CTCF Dnd41  *
#2 chr21 11619 11929 CTCF Dnd41  *
#3 chr21 11668 11978 CTCF Dnd41  *
#4 chr21 13183 13396 CTCF Dnd41  *
#5 chr21 14423 14580 CTCF Dnd41  *
#6 chr21 18733 18929 CTCF Dnd41  *

if(nrow(bed)>0){
bed$V2<-as.numeric(bed$V2)
bed$V3<-as.numeric(bed$V3)
bed<-unique(data.frame(chr=bed$V1,start=bed$V2,end=bed$V3,strand="*",antibody=bed$V4,celltype=bed$V5))
gr_range<-makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE) 
overlap<-countOverlaps(gr_cpg,gr_range)
print(length(which(overlap==1)))
w<-which(overlap>0)
if(length(w)>0){overlap[w]<-1}
}

if(nrow(bed)==0){overlap<-rep(0,nrow(r))}

ov<-data.frame(ov,overlap)

}
#names(ov)<-c("cpg",n.out)
names(ov)<-c("cpg",l2)

#cols <- names(ov)[-1]
#ov$ann <- apply( ov[ , cols ] , 1 , paste , collapse = "" )
#ov <- ov[ , !(names(ov)%in%cols) ]

#ov$ann <- do.call(paste, c(ov[cols], sep=""))

#ov2<-unite(ov, ann, -1,sep="") 
#a<-apply(ov2,1,function(x) nchar(x[2]))
m<-match(ov$cpg, df$cpg)
ov<-data.frame(df[m,],ov)
save(ov,file=paste0(chr,".Robj"))
#
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library(GenomicRanges)
library(data.table)

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/25states/25states_files.RData")
f<-ret
spl<-do.call("rbind",strsplit(f$filename,split="_"))
f$dataSource<-spl[,1]
f$seg_code<-gsub(".bed","",spl[,5])

state<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T)
state$STATE<-paste0("E",state$STATE)
m<-match(f$seg_code,state$STATE)
f$seg_explanation<-state[m,"NO."]

tiss<-read.table("jul2013.roadmapData_tissues.txt",sep="\t",he=T)
m<-match(f$dataSource,as.character(tiss$ID),)
f$tissue<-tiss[m,"Tissue"]
w<-which(f$tissue=="BLOOD")

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/25states/25states.RData")

#gr_range = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre)))
Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
Illumina450_dt=as.data.table(Illumina450)
Illumina450_dt[,cpgID:=row.names(Illumina450),]

gr_cpg = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

#overlap<-findOverlaps(gr_cpg,ret)
#overlap=as.data.table(findOverlaps(gr_cpg, ret))
#overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]
#Illumina450_dt[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
#Illumina450_dt[,cpgend:=end(gr_range[overlap_red$subjectHit])]
#Illumina450_dt[,NsubjectHits:=overlap_red$NsubjectHits]

ov<-as(findOverlaps(gr_cpg,ret), "List")
ov<-data.frame(ov)
ov<-ov[which(ov$value%in%w),]
df<-data.frame(ov,Illumina450_dt[ov$group])
df<-data.frame(df,f[df$value])
st<-unique(data.frame(cpg=df$cpgID,state=df$seg_explanation,tissue=df$tissue))
save(st,file="cpgbysegstate.Robj")





