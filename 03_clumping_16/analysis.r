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




clumped$rsq <- 2 * clumped$Effect^2 * clumped$Freq1 * (1 - clumped$Freq1)


cisthresh <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
transthresh <- c(1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14)
xc <- filter(clumped, cis)
xt <- filter(clumped, !cis)

# Count number of mQTLs

ar <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	ar[i] <- sum(xc$pval < cisthresh[i])
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- sum(xt$pval < transthresh[i])
}
mqtl_counts <- bind_rows(
	tibble(thresh=cisthresh, count=ar, cis='Cis'),
	tibble(thresh=transthresh, count=art, cis='Trans')
)

# Count number of mQTLs

arr <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	arr[i] <- sum(xc$PvalueARE < cisthresh[i])
}
arrt <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	arrt[i] <- sum(xt$PvalueARE < transthresh[i])
}
mqtl_countsr <- bind_rows(
	tibble(thresh=cisthresh, count=arr, cis='Cis'),
	tibble(thresh=transthresh, count=arrt, cis='Trans')
)

# Count number of mQTLs
arr <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	arr[i] <- sum(xc$PvalueMRE < cisthresh[i])
}
arrt <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	arrt[i] <- sum(xt$PvalueMRE < transthresh[i])
}
mqtl_countsrm <- bind_rows(
	tibble(thresh=cisthresh, count=arr, cis='Cis'),
	tibble(thresh=transthresh, count=arrt, cis='Trans')
)


##
y<-meffil.get.features("450k")
ncpg<-length(unique(y$name))


####
#n_independent_regions <- 1000000 used in Frank Dudbridge paper

n_independent_regions <- 1000000
#3 billion basepairs residing in 23 pairs of chromosomes
n_bases <- 3000000000

#SNP-CpG distance is 1 Mb 
cis_window <- 2000000

n_independent_regions_cis <- n_independent_regions / (n_bases / cis_window)
n_independent_regions_trans <- n_independent_regions - n_independent_regions_cis

#number of analysed CpGs - all probe on 450k array
ncpg <- retaincpg

ntest_cis <- ncpg * n_independent_regions_cis
ntest_trans <- ncpg * n_independent_regions_trans

# we used pval threshold of 1e-05
exp_cis <- ntest_cis * 1e-4
exp_trans <- ntest_trans * 1e-14
exp_cis
exp_trans

###
cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"
ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda.txt")
names(ss)<-c("study","nsamples16","lambda16")
ss[,1]<-gsub("_16","",ss[,1])
ss[,1]<-gsub("00_ARIES","ARIES",ss[,1])
#w<-which(ss[,1]%in%c("DunedinAge38"))
#if(length(w)>0){
#ss<-ss[-w,]}
ss4<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda4.txt")
names(ss4)<-c("study","nsamples04","lambda04")
ss4[,1]<-gsub("_04","",ss4[,1])
ss4[,1]<-gsub("00_ARIES","ARIES",ss4[,1])
w<-which(ss[,1]%in%c("ccg_old"))
if(length(w)>0){
ss<-ss[-w,]}

ss<-rbind(ss,c("SYS","NA","NA"),c("DunedinAge38","NA","NA"))
o<-order(ss[,1])
ss<-ss[o,]

ss$study_paper<-ss$study
ss$study_paper<-gsub("InterAct","EPIC_Norfolk",ss$study_paper)
ss$study_paper<-gsub("MARS_omni","MARS",ss$study_paper)
ss$study_paper<-gsub("Project_MinE_s27","MinE",ss$study_paper)
ss$study_paper<-gsub("as_cc","EGC_asthma",ss$study_paper)
ss$study_paper<-gsub("ccg","EGC_CTG",ss$study_paper)



w<-which(ss4[,1]%in%c("ccg_old","UCL_MRC_SCZ"))
if(length(w)>0){
ss4<-ss4[-w,]}


###
ss$nsamples01<-NA
ss$nsnp<-NA
ss$ncpg<-NA
ss$covariates<-NA
ss$males<-NA
ss$age<-NA

for (i in 1:nrow(ss)){
cat(ss[i,1],"\n")
load(paste0(cohort_dir,ss[i,1],"_01/results/01/cohort_descriptives.RData"))
ss$nsamples01[i]<-cohort_summary$mqtl_sample_size
ss$nsnp[i]<-cohort_summary$n_snp
ss$ncpg[i]<-cohort_summary$n_CpGs
ss$covariates[i]<-paste(cohort_summary$covariates,collapse=",")
ss$covariates[i]<-gsub("_numeric","",ss$covariates[i])
ss$covariates[i]<-gsub("_factor","",ss$covariates[i])
ss$cellcounts[i]<-paste(cohort_summary$predicted_cellcounts,collapse=",")
ss$predicted_cellcounts_type[i]<-cohort_summary$predicted_cellcounts_type
ss$males[i]<-round(cohort_summary$mqtl_n_males/cohort_summary$mqtl_sample_size,2)
ss$age[i]<-paste(round(cohort_summary$mqtl_mean_age,2)," (",round(cohort_summary$mqtl_min_age,2),"-",round(cohort_summary$mqtl_max_age,2),")",sep="")

if(cohort_summary$predicted_cellcounts_type=="NULL" & cohort_summary$predicted_cellcounts=="NULL")
{
ss$cellcounts[i]<-paste("Bcell,CD4T,CD8T,Eos,Mono,Neu,NK",sep="")
ss$predicted_cellcounts_type[i]<-"houseman"

}
}

m<-match(ss[,1],ss4[,1])
ss<-data.frame(ss,ss4[m,])

write.table(ss,"/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/descriptives.phase2.txt",sep="\t",quote=F,row.names=F,col.names=T)
#check samplesize differences between script 01 and 16
diff<-(ss$n-ss$nsamples01)
which(diff>0)

#total number of samples
sum(ss$n)
#[1] 28073


p1<-ggplot(ss,aes(x=nsamples01,y=nsnp/1000000))+
geom_point(aes(colour=factor(ss$study_paper),size=ss$nsamples01)) +
labs(x="Cohort N",y="N SNPs*1000000",size="N",colour="Study") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(plot=p1, file="../images/SNPsbyNcohort.pdf", width=7, height=7)

p2<-ggplot(ss,aes(x=nsamples01,y=ncpg/100000))+
geom_point(aes(colour=factor(ss$study_paper),size=ss$nsamples01)) +
labs(x="Cohort N",y="N CpGs*100000",size="N",colour="Study") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(plot=p2, file="../images/CpGsbyNcohort.pdf", width=7, height=7)

##
p3<-ggplot(ss,aes(x=nsamples16,y=lambda16))+
geom_point(aes(colour=factor(ss$study_paper),size=ss$nsamples16)) +
labs(x="Cohort N",y="lambda",size="N",colour="Study") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(plot=p3, file="../images/lambda_cisadj_byNcohort.pdf", width=7, height=7)


p4<-ggplot(ss,aes(x=nsamples04,y=lambda04))+
geom_point(aes(colour=factor(ss$study_paper),size=ss$nsamples04)) +
labs(x="Cohort N",y="lambda",size="N",colour="Study") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(plot=p4, file="../images/lambdabyNcohort.pdf", width=7, height=7)

p5<-ggplot(clumped,aes(x=(clumped$HetDf+1)))+
geom_histogram(binwidth=1)
ggsave(plot=p5, file="../images/numberofstudiesbymqtl.pdf", width=7, height=7)




##



y<-meffil.get.features("450k")
y<-y[which(!is.na(y$chromosome)),]
w<-which(y$chromosome=="chrY"|y$chromosome=="chrX")
cpgs<-unique(y$name)
sd.out<-data.frame(cpgs)
mean.out<-data.frame(cpgs)
for (i in 1:nrow(ss)){
cat(ss[i,1],"\n")
load(paste0(cohort_dir,ss[i,1],"_01/results/01/methylation_summary.RData"))
m<-match(cpgs,row.names(meth_summary))
sd<-meth_summary[m,"sd"]
sd.out<-data.frame(sd.out,sd)
mean.df<-meth_summary[m,"mean"]
mean.out<-data.frame(mean.out,mean.df)

}
row.names(sd.out)<-cpgs
row.names(mean.out)<-cpgs

sd.out<-sd.out[,-1]
names(sd.out)<-ss[,1]

mean.out<-mean.out[,-1]
names(mean.out)<-ss[,1]


for (i in 1:ncol(sd.out)){
cat(ss[i,1],"\n")
cat(length(which(is.na(sd.out[,i]))),"\n")
}

for (i in 1:ncol(sd.out)){
cat(ss[i,1],"\n")
cat(length(which(is.na(sd.out[-w,i]))),"\n")
}



sdmean<-data.frame(cpgs,sd=rowMeans(sd.out, na.rm=TRUE))
meanmean<-data.frame(cpgs,sd=rowMeans(mean.out, na.rm=TRUE))

probe.rm<-which(y$chromosome=="chrY"|y$chromosome=="chrX")
probe.rm<-unique(y[probe.rm,"name"])
w<-which(row.names(sd.out)%in%probe.rm)

studysd<-data.frame(study=ss[,1],sd=colMeans(sd.out[-w,], na.rm=TRUE))
studymean<-data.frame(study=ss[,1],mean=colMeans(mean.out[-w,], na.rm=TRUE))
write.table(studysd,"/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/sd.probes.txt",sep=" ",quote=F,col.names=T,row.names=F)
write.table(studymean,"/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/mean.probes.txt",sep=" ",quote=F,col.names=T,row.names=F)

#m<-match(row.names(sd.out),y$name)
#sd.out<-data.frame(y[m,c("chromosome","position")],sd.out)

m<-match(clumped$cpg,sdmean$cpgs)
clumped<-data.frame(clumped,cpgsd=sdmean[m,2])
m<-match(clumped$cpg,meanmean$cpgs)
clumped<-data.frame(clumped,cpgmean=meanmean[m,2])

sd.out.m<-melt(sd.out)
names(sd.out.m)<-c("study","sd")
p1 <- ggplot(sd.out.m, aes(x=as.factor(study), y=sd)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_text(size=6)) +
labs(x="study", y="sd")
ggsave(plot=p1, file="../images/sdbycohort.pdf", width=7, height=7)

mean.out.m<-melt(mean.out)
names(mean.out.m)<-c("study","mean")
p1 <- ggplot(mean.out.m, aes(x=as.factor(study), y=mean)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_text(size=6)) +
labs(x="study", y="mean")
ggsave(plot=p1, file="../images/meanbycohort.pdf", width=7, height=7)


clumped$id<-paste(clumped$snp,clumped$cpg,sep="_")
clumped$studycount<-0
###
path="/panfs/panasas01/shared-godmc/counts_2017/combined"
l<-list.files(path=path,pattern=".ge1.2.allcohorts.txt.gz")

for (i in 1:length(l)){
cat(i,"\n")
#r<-read.table(paste0("/panfs/panasas01/sscm/epzjlm/repo/goya/godmc/processed_data/phase2/assoclist/assoclist_",i,".gz"))
r<-read.table(paste0(path,"/",l[i]))
r$id<-paste(r$V2,r$V3,sep="_")
w<-which(r$id%in%clumped$id)
r<-r[w,]
w<-which(clumped$id%in%r$id)
clumped$studycount[w]<-r$V1
}
save(clumped,file="../results/16/16_clumpedwithstudycount.rdata")

#which association fall out if you compare fixed vs random
length(which(clumped$pval<clumped$PvalueARE))
#314495
length(which(clumped$pval<clumped$PvalueMRE))
#315096
length(which(clumped$pval>clumped$PvalueARE))
#0
length(which(clumped$pval>clumped$PvalueMRE))
#0
length(which(clumped$PvalueARE<clumped$PvalueMRE))
#312493
length(which(clumped$PvalueARE>clumped$PvalueMRE))
#2596

clumped2<-clumped[which(clumped$pval<1e-14),]
length(which(clumped2$pval<clumped2$PvalueARE))
#267763
length(which(clumped2$pval<clumped2$PvalueMRE))
#268364
length(which(clumped2$pval>clumped2$PvalueARE))
#0
length(which(clumped2$pval>clumped2$PvalueMRE))
#0
length(which(clumped2$PvalueARE<clumped2$PvalueMRE))
#268108
length(which(clumped2$PvalueARE>clumped2$PvalueMRE))
#256
#w<-which(clumped2$PvalueARE>clumped2$PvalueMRE)
#clumped2$mre<-"FE + MRE"
#clumped2$mre[w]<-"MRE only"
w<-which(clumped2$pval==0)
clumped2$pval[w]<-rep(1e-323,length(w))
w<-which(clumped2$PvalueARE==0)
clumped2$PvalueARE[w]<-rep(1e-323,length(w))
w<-which(clumped2$PvalueMRE==0)
clumped2$PvalueMRE[w]<-rep(1e-323,length(w))
clumped2$pval<-as.numeric(clumped2$pval)
clumped2$PvalueARE<-as.numeric(clumped2$PvalueARE)
clumped2$PvalueMRE<-as.numeric(clumped2$PvalueMRE)

w<-which(clumped2$pval<1e-14 & clumped2$PvalueMRE>1e-14)
clumped2$mre<-"FE + MRE +"
clumped2$mre[w]<-"FE + MRE - "

mean(abs(clumped2$Effect))
#[1] -0.009308104
mean(clumped2$HetDf)
#[1] 18.06114
mean(clumped2$StdErr)
#[1] 0.0171036
mean(clumped2$HetISq)
#[1] 52.77369
table(clumped2$HetDf)
#    0     1     2     3     4     5     6     7     8     9    10    11    12 
#  712   785  1167  1432  1803  1465  1798  2581  3629  4048  4779  4182  5206 
#   13    14    15    16    17    18    19    20    21    22    23 
# 6806  9311 12387 15769 20539 26034 28735 29971 32496 41776 31386 
mean(clumped2$tausq)
#[1] 0.01432297
mean(clumped2$TotalSampleSize)
#[1] 16615.94
table(clumped2$snptype)
#INDEL    SNP 
# 11292 277505 
table(clumped2$cis)
#FALSE   TRUE 
# 35966 252831 
####mean(clumped2$rsq)
####[1] 0.0425435
mean(clumped2$maf)
#[1] 0.2382847
mean(clumped2$cpgsd)
#[1] 0.04757534
mean(as.numeric(clumped2$studycount))
#[1] 6.1695

levels(clumped2$mre) <- rev(levels(clumped2$mre<-factor(clumped2$mre)))

p1 <- ggplot(clumped2, aes(x=as.factor(mre), y=-log10(pval))) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="-log10 FE pvalue")

p2 <- ggplot(clumped2, aes(x=as.factor(mre), y=-log10(PvalueMRE))) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="-log10 MRE pvalue")

p3 <- ggplot(clumped2, aes(x=as.factor(mre), y=abs(Effect))) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="abs FE Effect")

p4<- ggplot(clumped2, aes(x=as.factor(mre), y=StdErr)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="Std Err")

clumped2$nstudies<-clumped2$HetDf+1
p5 <- ggplot(clumped2, aes(x=as.factor(mre), y=nstudies)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="Number of phase2 Studies")

clumped2$studycount<-as.numeric(clumped2$studycount)
p6 <- ggplot(clumped2, aes(x=as.factor(mre), y=studycount)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="Number of phase1 Studies")

p7 <- ggplot(clumped2, aes(x=as.factor(mre), y=HetISq)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="HetISq")

p8 <- ggplot(clumped2, aes(x=as.factor(mre), y=tausq)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="tausq")

p9 <- ggplot(clumped2, aes(x=as.factor(mre), y=TotalSampleSize)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="TotalSampleSize")

p10 <- ggplot(clumped2, aes(x=as.factor(mre), y=rsq)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="rsq")

#p11 <- ggplot(clumped2, aes(x=as.factor(mre), y=maf)) +
#geom_boxplot() +
#theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
#labs(x="mre", y="maf")

p11 <- ggplot(clumped2, aes(x=as.factor(mre), y=cpgmean)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="cpgmean")

p12 <- ggplot(clumped2, aes(x=as.factor(mre), y=cpgsd)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(size=6),axis.text.y=element_text(size=6)) +
labs(x="mre", y="cpgsd")

pdf("../images/mrevsfe.pdf", width=7, height=7)
#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3,nrow=4)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3,nrow=4)
dev.off()




###
p1 <- ggplot(mqtl_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 38 cohorts")
ggsave(plot=p1, file="../images/mqtl_counts_thresholds_FE.pdf", width=7, height=7)

p1 <- ggplot(mqtl_countsr, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 38 cohorts")
ggsave(plot=p1, file="../images/mqtl_counts_thresholds_ARE.pdf", width=7, height=7)

p1 <- ggplot(mqtl_countsrm, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 38 cohorts")
ggsave(plot=p1, file="../images/mqtl_counts_thresholds_MRE.pdf", width=7, height=7)

###
temp4 <- subset(clumped, pval < 1e-14)
#b <- str_count(temp4$Direction, "\\?")
temp4$studycount<-temp4$HetDf+1
p1<-ggplot(temp4, aes(x=studycount)) + 
geom_histogram(binwidth=1)
ggsave(plot=p1, file="../images/nstudiesbydirection.pdf", width=7, height=7)

dim(temp4)
#[1] 288797     29
length(which(temp4$studycount<10))
#[1] 15372
length(which(temp4$studycount<5))
#[1] 4096
length(which(temp4$TotalSampleSize<10000))
#[1] 22159
length(which(temp4$TotalSampleSize<5000))
#[1] 6024
length(which(temp4$studycount<5&temp4$TotalSampleSize<5000))
#[1] 3989
length(which(temp4$studycount<5&temp4$TotalSampleSize<10000))
#[1] 4096

## Overlaps between fixed and random

arf <- rep(0, length(cisthresh))
arr <- rep(0, length(cisthresh))
arb <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	arf[i] <- sum(xc$pval < cisthresh[i] & xc$PvalueRandom > cisthresh[i])
	arr[i] <- sum(xc$pval > cisthresh[i] & xc$PvalueRandom < cisthresh[i])
	arb[i] <- sum(xc$pval < cisthresh[i] & xc$PvalueRandom < cisthresh[i])
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- sum(xt$pval < transthresh[i])
}

artf <- rep(0, length(transthresh))
artr <- rep(0, length(transthresh))
artb <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	artf[i] <- sum(xt$pval < transthresh[i] & xt$PvalueRandom > transthresh[i])
	artr[i] <- sum(xt$pval > transthresh[i] & xt$PvalueRandom < transthresh[i])
	artb[i] <- sum(xt$pval < transthresh[i] & xt$PvalueRandom < transthresh[i])
}


## Count number of CpGs

ar <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	ar[i] <- length(unique(subset(xc, pval < cisthresh[i])$cpg))
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- length(unique(subset(xt, pval < transthresh[i])$cpg))
}
cpg_counts <- bind_rows(
	tibble(thresh=cisthresh, count=ar, cis='Cis'),
	tibble(thresh=transthresh, count=art, cis='Trans')
)
p1 <- ggplot(cpg_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="CpGs with at least one mQTL from meta analysis of 7 cohorts")
ggsave(plot=p1, file="../images/cpg_counts_thresholds.pdf", width=7, height=7)


## Number of independent SNPs per cis and trans

sig <- subset(clumped, (cis & pval < 1e-7) | (!cis & pval < 1e-14))
clump_counts <- dplyr::group_by(sig, cpg, cis) %>%
	dplyr::summarise(n=n())

p1 <- ggplot(clump_counts, aes(x=as.factor(n))) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per CpG")
ggsave(plot=p1, file="../images/clump_counts_cpg.pdf", width=7, height=7)


## Number of hits per SNP

clump_counts_snp <- dplyr::group_by(sig, snp, cis) %>%
	dplyr::summarise(n=n()) %>%
	dplyr::group_by(n, cis) %>%
	dplyr::summarise(count=n())

p1 <- ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per SNP")
ggsave(plot=p1, file="../images/clump_counts_snp.pdf", width=10, height=7)



## Estimate Rsq

# What proportion of rsq is cis vs trans

clumped$rsq <- 2 * clumped$Effect^2 * clumped$Freq1 * (1 - clumped$Freq1)

rsq <- filter(clumped, pval < 5e-8) %>%
	dplyr::group_by(cpg, cis) %>%
	dplyr::summarise(rsq=sum(rsq))

rsqt <- filter(clumped, pval < 5e-8) %>%
	dplyr::group_by(cpg) %>%
	dplyr::summarise(rsq=sum(rsq))
rsq

p1 <- ggplot(rsq, aes(x=rsq)) +
geom_density(aes(fill=cis), alpha=0.3) +
labs(x="Rsq ")
ggsave(plot=p1, file="../images/rsq_density_cistrans.pdf", width=12, height=8)

# ggplot(rsqt, aes(x=rsq)) +
# geom_density() +
# labs(x="Rsq ")



# Plot cis against trans
# temp <- spread(rsq, key=cis, value=rsq)
# names(temp) <- c("cpg", "cis", "trans")
# temp <- subset(temp, !is.na(cis) & !is.na(trans))
# temp$tot <- temp$cis + temp$trans
# temp$rat <- temp$cis / temp$tot
# ggplot(temp, aes(x=cis, y=trans)) +
# geom_point(alpha=0.06) +
# labs(x="Rsq in cis", y="Rsq in trans")
# ggsave("../images/rsqcis_vs_rsqtrans.png", width=7, height=7)

# ggplot(temp, aes(x=rat)) +
# geom_density() +
# labs(x="Rsq_cis / (Rsq_cis + Rsq_trans)")
# ggsave("../images/rsqcis_vs_rsqtrans_ratio.png", width=7, height=7)




# mean(rsqt$rsq)

# median(subset(rsqt, cpg %in% temp$cpg)$rsq)
# median(subset(rsqt, ! cpg %in% temp$cpg)$rsq)

# group_by(rsq, cis) %>% summarise(rsq=sum(rsq)/450000)

#median()

## Position of CpG relative to SNP

temp1 <- subset(clumped, cpgchr == snpchr)
temp1$posdif <- temp1$snppos - temp1$cpgpos 
cisdist<-temp1[which(temp1$cis=="TRUE"),]
mediandist<-round(median(abs(cisdist$posdif))/1000,0)

p1 <- ggplot(subset(temp1, abs(posdif) < 2000000), aes(x=posdif)) +
	geom_density() +
	labs(x="Distance of SNP from CpG") +
	annotate(geom="text", x=0, y=4e-5, label=paste("Median distance = ",mediandist,"kb",sep=""), color="black")
ggsave(plot=p1, file="../images/snp_cpg_distance.pdf")

## Plot maf

clumped$maf <- clumped$Freq1
clumped$maf[clumped$maf > 0.5] <- 1 - clumped$maf[clumped$maf > 0.5]
p1 <- ggplot(clumped, aes(x=maf, y=abs(Effect))) +
geom_point(alpha=0.04) +
scale_x_log10() +
scale_y_log10()
ggsave(plot=p1, file="../images/maf_vs_beta.png", width=7, height=7)


## Directions

dir_count <- dplyr::group_by(clumped, Direction) %>%
	dplyr::summarise(count=n()) %>%
	filter(count > 100)


p1 <- ggplot(dir_count, aes(x=Direction, y=count)) +
geom_bar(stat="identity") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
ggsave(plot=p1, file="../images/directions_all.pdf", width=8, height=6)

## hterogeneity by direction

# temp3 <- subset(clumped, Direction %in% dir_count$Direction)
# ggplot(temp3, aes(x=Direction, y=-log10(HetPVal))) +
# geom_boxplot(fill="red") +
# theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

# temp33 <- group_by(temp3, Direction) %>% summarise(p = sum(HetPVal < 0.01)/n())

# ggplot(temp33, aes(x=Direction, y=p)) + geom_point() +
# theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))



temp4 <- subset(clumped, pval < 1e-14 & Direction %in% dir_count$Direction)
p1 <- ggplot(temp4, aes(x=Direction, y=HetISq)) +
geom_boxplot(fill="red") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Isq for mQTL with p < 1e-14")
ggsave(plot=p1, file="../images/isq_direction.pdf", width=10, height=6)

p1 <- ggplot(temp4, aes(x=Direction, y=tausq)) +
geom_boxplot(fill="red") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Tausq for mQTL with p < 1e-14") + scale_y_log10()
ggsave(plot=p1, file="../images/tausq_direction.pdf", width=10, height=6)

p1 <- ggplot(temp4, aes(x=cpgsd, y=tausq)) +
geom_point(alpha=0.04) +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Tausq for mQTL with p < 1e-14") + scale_y_log10()
labs(x="SD for mQTL CpG")
ggsave(plot=p1, file="../images/tausq_cpgsd.pdf", width=10, height=6)

temp4$sdcat<-cut(as.numeric(as.character(temp4$cpgsd)), breaks = seq(0,1,by=0.01))
labs <- data.frame(table(temp4$sdcat))
m<-match(temp4$sdcat,labs[,1])
temp4<-data.frame(temp4,Ncatsd=as.character(labs[m,-1]))

m <- data.frame(c(by(abs(temp4$tausq), temp4$sdcat, max)))
m2<-match(temp4$sdcat,row.names(m))
temp4<-data.frame(temp4,maxtausq=abs(m[m2,]))

p1<-ggplot(temp4, aes(x=sdcat, y=tausq)) +
geom_boxplot() +
geom_text(data=temp4, aes(x=temp4$sdcat,y = (maxtausq+0.1),label = Ncatsd),vjust = 0,size=3) +
ylab("Tau^2") 
ggsave(p1,file="../images/sdcatvstausq.pdf",height=6,width=16)

p1 <- ggplot(temp4, aes(x=HetISq, y=tausq)) +
geom_point(alpha=0.04) +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Tausq for mQTL with p < 1e-14") + scale_y_log10()
labs(x="Isq for mQTL with p < 1e-14")
ggsave(plot=p1, file="../images/tausq_isq.pdf", width=10, height=6)

temp4$Effect_abs<-abs(temp4$Effect)
temp4$i2cat<-cut(as.numeric(as.character(temp4$HetISq)), breaks = seq(0,100,by=10))
w<-which(is.na(temp4$i2cat))
table(temp4$HetISq[w])
temp4$i2cat[w]<-"(0,10]"

labs <- data.frame(table(temp4$i2cat))
m<-match(temp4$i2cat,labs[,1])
temp4<-data.frame(temp4,Ncati2=as.character(labs[m,-1]))

m <- data.frame(c(by(abs(temp4$Effect), temp4$i2cat, max)))
m2<-match(temp4$i2cat,row.names(m))
temp4<-data.frame(temp4,maxbeta=abs(m[m2,]))


p1<-ggplot(temp4, aes(x=i2cat, y=Effect_abs)) +
geom_boxplot() +
geom_text(data=temp4, aes(x=temp4$i2cat,y = (maxbeta+0.1),label = Ncati2),vjust = 0,size=3) +
ylab("Beta coefficient") 
ggsave(p1,file="../images/I2catvsbeta.pdf",height=6,width=16)

#
temp4$tau2cat<-cut(as.numeric(as.character(temp4$tausq)), breaks = seq(0,2,by=0.01))
w<-which(is.na(temp4$tau2cat))
table(temp4$tau2cat[w])
temp4$tau2cat[w]<-"(0,0.01]"

labs <- data.frame(table(temp4$tau2cat))
m<-match(temp4$tau2cat,labs[,1])
temp4<-data.frame(temp4,Ncattau2=as.character(labs[m,-1]))

m <- data.frame(c(by(abs(temp4$Effect), temp4$tau2cat, max)))
m2<-match(temp4$tau2cat,row.names(m))
temp4<-data.frame(temp4,maxbeta_tau=abs(m[m2,]))

p1<-ggplot(temp4, aes(x=tau2cat, y=Effect_abs)) +
geom_boxplot() +
geom_text(data=temp4, aes(x=temp4$tau2cat,y = (maxbeta_tau+0.01),label = Ncattau2),vjust = 0,size=3) +
ylab("Beta coefficient") 
ggsave(p1,file="../images/tau2catvsbeta.pdf",height=6,width=16)

temp4$studycount<-factor(temp4$studycount,levels=temp4$studycount[order(as.numeric(as.character(temp4$studycount)))])
p1<-ggplot(temp4, aes(x=studycount, y=tausq)) +
geom_boxplot() +
#geom_text(data=temp4, aes(x=temp4$tausq,y = (maxbeta_tau+0.01),label = Ncattau2),vjust = 0,size=3) +
ylab("Tau^2") 
ggsave(p1,file="../images/tau2vsstudycount.pdf",height=6,width=16)

p1<-ggplot(temp4, aes(x=studycount, y=HetISq)) +
geom_boxplot() +
#geom_text(data=temp4, aes(x=temp4$tausq,y = (maxbeta_tau+0.01),label = Ncattau2),vjust = 0,size=3) +
ylab("I2") 
ggsave(p1,file="../images/i2vsstudycount.pdf",height=6,width=16)

#
temp4$qcat<-cut(as.numeric(as.character(temp4$HetChiSq)), breaks = seq(0,2000,by=10))
w<-which(is.na(temp4$qcat))
#table(temp4$qcat[w])
temp4$qcat[w]<-"(0,10]"

labs <- data.frame(table(temp4$qcat))
m<-match(temp4$qcat,labs[,1])
temp4<-data.frame(temp4,Ncatq=as.character(labs[m,-1]))

m <- data.frame(c(by(abs(temp4$Effect), temp4$qcat, max)))
m2<-match(temp4$qcat,row.names(m))
temp4<-data.frame(temp4,maxbeta_q=abs(m[m2,]))

p1<-ggplot(temp4, aes(x=qcat, y=Effect_abs)) +
geom_boxplot() +
geom_text(data=temp4, aes(x=temp4$qcat,y = (maxbeta_q+0.1),label = Ncatq),vjust = 0,size=3) +
ylab("Beta coefficient") 
ggsave(p1,file="../images/qcatvsbeta.pdf",height=6,width=16)



# Conditioanl 

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_conditional.rdata")

## Number of independent SNPs per cis and trans

sig <- subset(conditional, (cis & pval < 1e-7) | (!cis & pval < 1e-14))
clump_counts <- dplyr::group_by(sig, cpg, cis) %>%
	dplyr::summarise(n=n())

p1 <- ggplot(clump_counts, aes(x=n)) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from conditional analysis (p < 1e-7; 1e-14)", y="mQTLs per CpG")
ggsave(p1, file="../images/conditional_counts_cpg.pdf", width=7, height=7)


## Number of hits per SNP

clump_counts_snp <- dplyr::group_by(sig, snp, cis) %>%
	dplyr::summarise(n=n()) %>%
	dplyr::group_by(n, cis) %>%
	dplyr::summarise(count=n())

p2 <- ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from conditional analysis (p < 1e-7; 1e-14)", y="mQTLs per SNP")
ggsave(p2, file="../images/conditional_counts_snp.pdf", width=7, height=7)




