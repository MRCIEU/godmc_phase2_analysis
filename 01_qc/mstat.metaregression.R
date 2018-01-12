library(ggplot2)
library(metafor)
library(dplyr)   # for sorting data.frames
library(getmstatistic)  # for calculating M statistics
library(gridExtra)

cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"

load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/mstat/mstats_chr20.RData")
dframe2<-unique(data.frame(study_names_in=dframe$study_names_in,M=dframe$M,M_sd=dframe$M_sd,M_se=dframe$M_se,dframe$study))
dframe2[,1]<-gsub("00_ARIES","ARIES",dframe2[,1])
path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis"
w<-which(dframe2[,1]%in%c("DunedinAge38"))
if(length(w)>0){
dframe2<-dframe2[-w,]}

dframe2$age<-NA
dframe2$nsnps<-NA
dframe2$relatedness<-NA
dframe2$ncpg<-NA
dframe2$nmales<-NA
dframe2$maf<-NA
dframe2$info<-NA
dframe2$ccmethod<-NA
dframe2$bmi<-NA
dframe2$n_bmi<-NA
dframe2$height<-NA
dframe2$sd_bmi<-NA
dframe2$sd_height<-NA


for (i in 1:nrow(dframe2)){
cat(dframe2[i,1],"\n")
load(paste0(cohort_dir,dframe2[i,1],"_01/results/01/cohort_descriptives.RData"))
age<-cohort_summary$mqtl_mean_age
dframe2$age[i]<-age
nsnps<-cohort_summary$n_snp
dframe2$nsnps[i]<-nsnps
ncpgs<-cohort_summary$n_CpGs
dframe2$ncpg[i]<-ncpgs
nmales<-cohort_summary$mqtl_n_males/cohort_summary$mqtl_sample_size
dframe2$nmales[i]<-nmales
maf<-cohort_summary$mean_maf
dframe2$maf[i]<-maf
info<-cohort_summary$mean_info
dframe2$info[i]<-info
ccmethod<-cohort_summary$predicted_cellcounts_type
dframe2$ccmethod[i]<-ccmethod


if(length(which(names(cohort_summary)%in%c("mean_BMI")))==1){
bmi<-cohort_summary$mean_BMI
dframe2$bmi[i]<-bmi

sd_bmi<-cohort_summary$sd_BMI
dframe2$sd_bmi[i]<-sd_bmi
}

if(length(which(names(cohort_summary)%in%c("BMI_sample_size")))==1){
n_bmi<-cohort_summary$BMI_sample_size
dframe2$n_bmi[i]<-n_bmi

sd_bmi<-cohort_summary$sd_BMI
dframe2$sd_bmi[i]<-sd_bmi
}

if(length(which(names(cohort_summary)%in%c("mean_Height")))==1){
height<-cohort_summary$mean_Height
dframe2$height[i]<-height
sd_height<-cohort_summary$sd_Height
dframe2$sd_height[i]<-sd_height

}

}

w<-which(dframe2$ccmethod=="NULL")
dframe2$ccmethod[w]<-"houseman"
w<-which(dframe2$ccmethod=="bakulski")
dframe2$ccmethod[w]<-"cord"
dframe2$ccmethod<-as.factor(dframe2$ccmethod)
dframe2$ccmethod = with(dframe2, factor(dframe2$ccmethod, levels = rev(levels(dframe2$ccmethod))))

w<-which(dframe2$age==0)
dframe2$height[w]<-NA
dframe2$sd_height[w]<-NA
dframe2$bmi[w]<-NA
dframe2$sd_bmi[w]<-NA

ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda4.txt")
names(ss)<-c("study","n","lambda4")
ss[,1]<-gsub("_04","",ss[,1])
ss[,1]<-gsub("00_ARIES","ARIES",ss[,1])
m<-match(dframe2$study_names_in,ss[,1])
dframe2$samplesize<-ss[m,"n"]
dframe2$lambda4<-ss[m,"lambda4"]

ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda.txt")
names(ss)<-c("study","n","lambda16")
ss[,1]<-gsub("_16","",ss[,1])
ss[,1]<-gsub("00_ARIES","ARIES",ss[,1])
m<-match(dframe2$study_names_in,ss[,1])
dframe2$samplesize<-ss[m,"n"]
dframe2$lambda16<-ss[m,"lambda16"]

w<-which(dframe2$study_names_in%in%c("BorninBradford","Project_MinE_s27","GLAKU"))
dframe2$chip<-"450k"
dframe2$chip[w]<-"EPIC"

rel<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt")
rel[,1]<-gsub("_01","",rel[,1])
m<-match(dframe2$study_names_in,rel[,1])
dframe2$relatedness<-rel[m,-1]
dframe2$relatedness<-factor(dframe2$relatedness)
dframe2$relatedness = with(dframe2, factor(dframe2$relatedness, levels = rev(levels(dframe2$relatedness))))

ancestry<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortancestry.txt",sep="\t")
ancestry[,1]<-gsub("MARS","MARS_omni",ancestry[,1])
m<-match(dframe2$study_names_in,ancestry[,1])
dframe2$ancestry<-ancestry[m,-1]
dframe2$ancestry<-factor(dframe2$ancestry)
dframe2$ancestry = with(dframe2, factor(dframe2$ancestry, levels = rev(levels(dframe2$ancestry))))

dframe2$ancestry_ns<-as.character(dframe2$ancestry)
w<-which(dframe2$ancestry_ns%in%c("UK","The Netherlands","Sweden","Scotland","Germany","Canada","Denmark","Estonia","Finland"))
dframe2$ancestry_ns[w]<-"north"
w<-which(dframe2$ancestry_ns%in%c("Australia","New Zealand","France","Spain","Italy"))
dframe2$ancestry_ns[w]<-"south"
dframe2$ancestry_ns<-as.factor(dframe2$ancestry_ns)

sd.p<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/sd.probes.txt",sep=" ",he=T)
m<-match(dframe2$study_names_in,sd.p[,1])
dframe2$sd.probe<-sd.p[m,-1]

mean.p<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/mean.probes.txt",sep=" ",he=T)
m<-match(dframe2$study_names_in,mean.p[,1])
dframe2$mean.probe<-mean.p[m,-1]

norm.method<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortnormalizationmethod.txt",sep="\t",he=T)
m<-match(dframe2$study_names_in,norm.method[,1])
dframe2$norm.method<-norm.method[m,-1]
w<-which(dframe2$norm.method%in%c("BMIQ","CPACOR","Dasen","glm","methylumi","SWAN"))
dframe2$norm.method_fn<-as.character(dframe2$norm.method)
dframe2$norm.method_fn[w]<-"Other"
dframe2$norm.method_fn<-as.factor(dframe2$norm.method_fn)

dframe2$chip<-as.factor(dframe2$chip)
dframe2$nsnps<-as.numeric(dframe2$nsnps)
dframe2$ncpg<-as.numeric(dframe2$ncpg)
dframe2$samplesize<-as.numeric(dframe2$samplesize)

#bmi missingness
miss<-data.frame(study=dframe2$study_names_in,n_bmi=dframe2$n_bmi,samplesize=dframe2$samplesize,perc=dframe2$n_bmi/dframe2$samplesize)
w<-which(miss$perc<0.85 & dframe2$study!="BorninBradford") #BiB are mums + cord, for cords there is no BMI

#PRECISESADS/ALS1/ALS2 don't have BMI or Height

if(length(w)>0){
dframe2$height[w]<-NA
dframe2$sd_height[w]<-NA
dframe2$bmi[w]<-NA
dframe2$sd_bmi[w]<-NA
}

dframe2$highbmi<-NA
w<-which(dframe2$bmi>25)
dframe2$highbmi[w]<-1
w<-which(dframe2$bmi<25)
dframe2$highbmi[w]<-0

dframe2$vi<-0
#res<-rma(yi=M,vi=M_sd^2,mods=~age+nsnps+samplesize+relatedness+ancestry,data=dframe2)
#res<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe+ancestry,data=dframe2)
#yi: vector of length k with the observed effect sizes or outcomes. See ‘Details’.
#vi: vector of length k with the corresponding sampling variances. See ‘Details’.
#res2<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+height+bmi+ccmethod+mean.probe+sd.probe,data=dframe2)
#res3<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe,data=dframe2)

print(rma(yi=M,vi=vi,mods=~age,data=dframe2))
print(rma(yi=M,vi=vi,mods=~nmales,data=dframe2))
print(rma(yi=M,vi=vi,mods=~ancestry_ns,data=dframe2))
print(rma(yi=M,vi=vi,mods=~nsnps,data=dframe2)) #0.0092
print(rma(yi=M,vi=vi,mods=~chip,data=dframe2))
print(rma(yi=M,vi=vi,mods=~ncpg,data=dframe2))
print(rma(yi=M,vi=vi,mods=~samplesize,data=dframe2))
print(rma(yi=M,vi=vi,mods=~relatedness,data=dframe2))
print(rma(yi=M,vi=vi,mods=~maf,data=dframe2))#0.0625
print(rma(yi=M,vi=vi,mods=~info,data=dframe2))
print(rma(yi=M,vi=vi,mods=~lambda4,data=dframe2))
print(rma(yi=M,vi=vi,mods=~lambda16,data=dframe2))
print(rma(yi=M,vi=vi,mods=~ccmethod,data=dframe2))
print(rma(yi=M,vi=vi,mods=~mean.probe,data=dframe2))
print(rma(yi=M,vi=vi,mods=~sd.probe,data=dframe2))#<.0001
print(rma(yi=M,vi=vi,mods=~ancestry,data=dframe2))
print(rma(yi=M,vi=vi,mods=~height,data=dframe2))#
print(rma(yi=M,vi=vi,mods=~bmi,data=dframe2))#0.0283
print(rma(yi=M,vi=vi,mods=~norm.method,data=dframe2))
print(rma(yi=M,vi=vi,mods=~sd_bmi,data=dframe2))#0.6102
print(rma(yi=M,vi=vi,mods=~sd_height,data=dframe2))
print(rma(yi=M,vi=vi,mods=~nsnps+bmi,data=dframe2))

#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub    
#intrcpt    2.0731  0.6883   3.0120  0.0026   0.7241  3.4221  **
#nsnps     -0.0000  0.0000  -1.8515  0.0641  -0.0000  0.0000   .
#bmi       -0.0458  0.0235  -1.9488  0.0513  -0.0918  0.0003   .

print(rma(yi=M,vi=vi,mods=~nsnps+bmi+sd_bmi,data=dframe2))
#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub     
#intrcpt    2.7020  0.6688   4.0402  <.0001   1.3912   4.0127  ***
#nsnps     -0.0000  0.0000  -2.6631  0.0077  -0.0000  -0.0000   **
#bmi       -0.0965  0.0298  -3.2399  0.0012  -0.1549  -0.0381   **
#sd_bmi     0.2261  0.0941   2.4017  0.0163   0.0416   0.4106    *

age18<-which(dframe2$age<18)
print(rma(yi=M,vi=vi,mods=~nsnps+bmi+sd_bmi,data=dframe2[-age18,]))
#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub    
#intrcpt    2.7321  0.9928   2.7518  0.0059   0.7861   4.6780  **
#nsnps     -0.0000  0.0000  -2.3835  0.0172  -0.0000  -0.0000   *
#bmi       -0.0956  0.0391  -2.4465  0.0144  -0.1721  -0.0190   *
#sd_bmi     0.2035  0.1113   1.8290  0.0674  -0.0146   0.4216   .

w<-which(dframe2$study_names_in%in%c("BAMSE"))
print(rma(yi=M,vi=vi,mods=~nsnps+bmi+sd_bmi,data=dframe2[-w,]))
#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub     
#intrcpt    2.9265  0.7915   3.6972  0.0002   1.3751   4.4779  ***
#nsnps     -0.0000  0.0000  -2.5493  0.0108  -0.0000  -0.0000    *
#bmi       -0.1042  0.0334  -3.1234  0.0018  -0.1696  -0.0388   **
#sd_bmi     0.2146  0.0982   2.1861  0.0288   0.0222   0.4070    *

#test for interaction
#split in high and normal BMI
#sd BMI against Mvalue

print(rma(yi=M,vi=vi,mods=~nsnps+bmi*sd_bmi,data=dframe2))

#Model Results:

#            estimate      se     zval    pval    ci.lb    ci.ub   
#intrcpt       2.5498  1.7293   1.4745  0.1404  -0.8396   5.9393   
#nsnps        -0.0000  0.0000  -2.5112  0.0120  -0.0000  -0.0000  *
#bmi          -0.0903  0.0713  -1.2675  0.2050  -0.2300   0.0493   
#sd_bmi        0.2687  0.4545   0.5911  0.5544  -0.6222   1.1596   
#bmi:sd_bmi   -0.0018  0.0187  -0.0959  0.9236  -0.0384   0.0348   


w<-which(dframe2$highbmi==1) #12
print(rma(yi=M,vi=vi,mods=~nsnps+sd_bmi,data=dframe2[w,]))

#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub    
#intrcpt    0.3133  0.4483   0.6989  0.4846  -0.5654   1.1921    
#nsnps     -0.0000  0.0000  -2.7583  0.0058  -0.0000  -0.0000  **
#sd_bmi     0.2177  0.1114   1.9544  0.0506  -0.0006   0.4360   .

w<-which(dframe2$highbmi==0) #10
print(rma(yi=M,vi=vi,mods=~nsnps+sd_bmi,data=dframe2[w,]))

#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub   
#intrcpt    2.7441  1.1158   2.4594  0.0139   0.5572   4.9310  *
#nsnps     -0.0000  0.0000  -2.4452  0.0145  -0.0000  -0.0000  *
#sd_bmi     0.0690  0.1127   0.6116  0.5408  -0.1520   0.2899   


print(rma(yi=M,vi=vi,mods=~nsnps+age+bmi+sd_bmi+sd.probe,data=dframe2))
#          estimate      se     zval    pval     ci.lb     ci.ub     
#intrcpt     2.8978  0.7946   3.6470  0.0003    1.3405    4.4552  ***
#nsnps      -0.0000  0.0000  -2.8513  0.0044   -0.0000   -0.0000   **
#age         0.0052  0.0049   1.0697  0.2848   -0.0043    0.0147     
#bmi        -0.0832  0.0416  -1.9970  0.0458   -0.1648   -0.0015    *
#sd_bmi      0.2035  0.0779   2.6107  0.0090    0.0507    0.3562   **
#sd.probe  -27.0831  7.2771  -3.7217  0.0002  -41.3459  -12.8203  ***

print(rma(yi=M,vi=vi,mods=~nsnps+age+bmi+sd_bmi+sd.probe,data=dframe2))
#Model Results:

#          estimate      se     zval    pval     ci.lb     ci.ub     
#intrcpt     2.8978  0.7946   3.6470  0.0003    1.3405    4.4552  ***
#nsnps      -0.0000  0.0000  -2.8513  0.0044   -0.0000   -0.0000   **
#age         0.0052  0.0049   1.0697  0.2848   -0.0043    0.0147     
#bmi        -0.0832  0.0416  -1.9970  0.0458   -0.1648   -0.0015    *
#sd_bmi      0.2035  0.0779   2.6107  0.0090    0.0507    0.3562   **
#sd.probe  -27.0831  7.2771  -3.7217  0.0002  -41.3459  -12.8203  ***

w<-which(is.na(dframe2$bmi))
dframe3<-dframe2[-w,]


cor.test(dframe3$sd_bmi,dframe3$sd.probe)
#p-value = 0.1745 cor=0.3003152 

cor.test(dframe3$sd_bmi,dframe3$bmi)
#      cor 
#0.7184677 
#t = 4.6281, df = 20, p-value = 0.0001657
rma(yi=M,vi=vi,mods=~nsnps+bmi+ancestry,data=dframe2)

rma(yi=bmi,vi=vi,mods=~sd_bmi+sd.probe,data=dframe3)


p1<-ggplot(dframe3,aes(x=bmi,y=sd_bmi))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average BMI") +
ylab("BMI sd") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/sdbmibyBMI.pdf",width=8,height=6)


p1<-ggplot(dframe3,aes(x=sd.probe,y=sd_bmi))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average sd probe") +
ylab("BMI sd") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/sdprobebyBMIsd.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=bmi,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average BMI") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/MvaluebyBMI.pdf",width=8,height=6)

p1<-ggplot(dframe2,aes(x=maf,y=M))+
geom_point(aes(colour=factor(dframe2$study_names_in),size=dframe2$samplesize)) +
xlab("average MAF") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/Mvaluebymaf.pdf",width=8,height=6)

p1<-ggplot(dframe2,aes(x=nsnps,y=M))+
geom_point(aes(colour=factor(dframe2$study_names_in),size=dframe2$samplesize)) +
xlab("nSNPs") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/Mvaluebynsnps.pdf",width=8,height=6)

p1<-ggplot(dframe2,aes(x=sd.probe,y=M))+
geom_point(aes(colour=factor(dframe2$study_names_in),size=dframe2$samplesize)) +
xlab("average SD across all probes") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/Mvaluebysd.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=sd_bmi,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("BMI SD") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="./images/Mvaluebysdbmi.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=sd_bmi,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average BMI SD") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold")) +
facet_grid(.~highbmi)
ggsave(p1,file="./images/Mvaluebysdbmi_highvslowBMI.pdf",width=8,height=6)

# Sort getmstatistic_results dataframe by M statistics
dframe$study_names_in<-gsub("00_ARIES","ARIES",dframe$study_names_in)
getm_res_srtd <- dplyr::arrange(dframe, M)

# First, drop duplicate study_names in the sorted getmsatistic_results
# choose a snp with the maximum number of studies
getm_res_srtd_nodups <- subset(getm_res_srtd, getm_res_srtd$snp == 6)

# Checking dimensions to confirm that we have 38 studies
#str(getm_res_srtd_nodups)

#
getm_res_plus_n <- merge(getm_res_srtd_nodups, ss, by.x="study_names_in", by.y="study")
#getm_res_plus_n <- merge(getm_res_srtd, ss, by.x="study_names_in", by.y="study")


# Sort data.frame by M statistics
getm_res_plus_n_srtd <- dplyr::arrange(getm_res_plus_n, M)
metafor_results_fe <- metafor::rma.uni(yi = getm_res_plus_n_srtd[, "M"], sei = getm_res_plus_n_srtd[, "M_se"], weighted = T, slab = getm_res_plus_n_srtd[, "study_names_in"], method = "FE")

# Compute inverse-variance weighted random effects model
metafor_results_dl <- metafor::rma.uni(yi = getm_res_plus_n_srtd[, "M"], sei = getm_res_plus_n_srtd[, "M_se"], weighted = T, knha = T, slab = getm_res_plus_n_srtd[, "study_names_in"], method = "DL")
metafor_results_dl


# Plotting:

# set margins
pdf("./images/mstat.forestplot.pdf",height=6,width=6)
par(mar=c(4,4,1,2))

#forest(metafor_results_fe,cex=0.01,psize=0.01)
# generate forest plot
forest(metafor_results_fe, xlim=c(-3, 1.6), at=c(-1, -0.5, 0, 0.5, 1), cex=0.66, xlab = "M statistic", ilab=getm_res_plus_n_srtd$n, ilab.xpos = c(-1.1), ilab.pos = c(2,2), addfit=F)

# Adding labels
text(1.6, length(getm_res_plus_n_srtd$study_names_in)+2, "M statistic [95% CI]", pos=2, cex=0.66)
text(-1.4, length(getm_res_plus_n_srtd$study_names_in)+2, "N", pos=4, cex=0.66)
#text(-1.39, 50, "Controls", pos=4, cex=0.66)
text(-3, length(getm_res_plus_n_srtd$study_names_in)+2, "Study", pos=4, cex=0.66)
dev.off()
