library(GenomicRanges)
library(data.table)
library(tidyr)
library(meffil)
library(SDMTools)
library(dplyr)

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

data$MAF<-data$Freq1
w<-which(data$MAF>0.5)
data$MAF[w]<-1-data$MAF[w]

w<-which(data$cpg_cis=="cis+trans" & data$cis!="TRUE")
data[w,"cpg_cis"]<-"cis_trans_trans"
w<-which(data$cpg_cis=="cis+trans" & data$cis=="TRUE")
data[w,"cpg_cis"]<-"cis_trans_cis"


data<-data.table(data)
data2<-data
data2$cpg_cis<-"All"
data2<-rbind(data,data2)

p1 <- ggplot(data2, aes(x=as.factor(cpg_cis), y=abs(MAF))) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="MAF")
ggsave(plot=p1, file="./images/assoc_MAF_cpg_ciscategory.pdf", width=7, height=7)


load("meancpg.Robj")
m<-match(data2$cpg,df2$cpg)
data2<-data.frame(data2,df2[m,c("meancpg","sdcpg")])

amb<-subset(data,cpg_cis=="cis_trans_trans"|cpg_cis=="cis_trans_cis" )

data$Direction<-""
w1<-which(data$Effect<0)
data$Direction[w1]<-"-"
w1<-which(data$Effect>0)
data$Direction[w1]<-"+"

a_data<-aggregate(Direction ~ cpg, data = data, paste, collapse = ",")
m<-match(a_data$cpg,data2$cpg)
a_data<-data.frame(a_data,data2[m,c("cpg_cis","meancpg","sdcpg")])
a_data$cpg_cis<-gsub("cis_trans_trans","cis_trans",a_data$cpg_cis)

table(a_data$Direction)

p1 <- ggplot(data2, aes(x=as.factor(cpg_cis), y=abs(Effect))) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="abs effect")
ggsave(plot=p1, file="./images/assoc_effect_cpg_ciscategory.pdf", width=7, height=7)


g1<-grep("+",a_data$Direction,fixed=T)
g2<-grep("-",a_data$Direction,fixed=T)
g<-intersect(g1,g2)
a_data$Direction2<-""
a_data$Direction2[g1]<-"+"
a_data$Direction2[g2]<-"-"
a_data$Direction2[g]<-"-+"
a_data%>%group_by(cpg_cis,Direction2)%>%summarize(mean=mean(meancpg))

a_data2<-a_data
a_data2$cpg_cis<-"All"
a_data2<-rbind(a_data,a_data2)

p1 <- ggplot(a_data2, aes(x=Direction2, y=meancpg)) +
geom_boxplot() +
facet_wrap(~cpg_cis,scales="free_y") +
labs(x="Direction", y="Weighted Mean methylation")
ggsave(plot=p1, file="./images/assoc_direction_meancpg_ciscategory.pdf", width=7, height=7)


p1 <- ggplot(data2, aes(x=as.factor(cpg_cis), y=abs(Effect))) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="abs effect")
ggsave(plot=p1, file="./images/assoc_effect_cpg_ciscategory.pdf", width=7, height=7)


p1 <- ggplot(data2, aes(x=as.factor(cpg_cis), y=HetISq)) +
geom_boxplot() +
theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="Category", y="I2")
ggsave(plot=p1, file="./images/assoc_i2_cpg_ciscategory.pdf", width=7, height=7)

data2%>%group_by(cpg_cis)%>%summarize(mean=mean(HetISq))
# A tibble: 5 x 2
#          cpg_cis     mean
#            <chr>    <dbl>
#1             All 45.26092
#2        cis only 45.95120
#3   cis_trans_cis 40.60503
#4 cis_trans_trans 35.84082
#5      trans only 52.73976

data2$cpg_cis2<-data2$cpg_cis
data2$cpg_cis2<-gsub("cis_trans_cis","cis_trans",data2$cpg_cis2)
data2$cpg_cis2<-gsub("cis_trans_trans","cis_trans",data2$cpg_cis2)

data2%>%group_by(cpg_cis2)%>%summarize(mean=mean(HetISq))
# A tibble: 4 x 2
#    cpg_cis2     mean
#       <chr>    <dbl>
#1        All 45.26092
#2   cis only 45.95120
#3  cis_trans 38.52523
#4 trans only 52.73976

data2%>%group_by(cpg_cis)%>%summarize(mean=mean(HetISq))
# A tibble: 5 x 2
#          cpg_cis     mean
#            <chr>    <dbl>
#1             All 45.26092
#2        cis only 45.95120
#3   cis_trans_cis 40.60503
#4 cis_trans_trans 35.84082
#5      trans only 52.73976

data2%>%group_by(cpg_cis)%>%summarize(mean=mean(MAF))
# A tibble: 5 x 2
#          cpg_cis      mean
#            <chr>     <dbl>
#1             All 0.2336860
#2        cis only 0.2331293
#3   cis_trans_cis 0.2217990
#4 cis_trans_trans 0.2446098
#5      trans only 0.2558984

data2$abs_Effect<-abs(data2$Effect)
# A tibble: 5 x 2
#          cpg_cis      mean
#            <chr>     <dbl>
#1             All 0.2586500
#2        cis only 0.2598534
#3   cis_trans_cis 0.2910761
#4 cis_trans_trans 0.1987415
#5      trans only 0.2568757



