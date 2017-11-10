#extract results files
#/panfs/panasas01/shared-godmc/sftp/GoDMC/ecarnero/TwinsUK/results/01
#cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
#results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"


#cd $results_dir

#for cohort in ${cohort_dir}*01.tgz
#do
#echo $cohort


#	cohortname=$(basename "$cohort" .tgz)
#	echo $cohortname
#	mkdir -p $cohort_dir/${cohortname}
#	cd $cohort_dir/${cohortname}


	#tar -xvf ../${cohortname}.tgz results/01/cohort_descriptives.RData
    #tar -xvf ../${cohortname}.tgz config
    #related=`grep "related=" config |perl -pe 's/related=//g'|perl -pe 's/"//g'`
    #echo $cohortname $related >>/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt 
 #   tar -xvf ../${cohortname}.tgz results/01/methylation_summary.RData
#done
###
library(ggplot2)
library(metafor)
library(dplyr)   # for sorting data.frames
library(getmstatistic)  # for calculating M statistics
library(gridExtra)

cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"

load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/mstat/mstats.RData")
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
dframe2$height<-NA


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
dframe2$bmi[i]<-bmi}

if(length(which(names(cohort_summary)%in%c("mean_Height")))==1){
height<-cohort_summary$mean_Height
dframe2$height[i]<-height
}

}


w<-which(dframe2$ccmethod=="NULL")
dframe2$ccmethod[w]<-"houseman"
w<-which(dframe2$ccmethod=="bakulski")
dframe2$ccmethod[w]<-"cord"
dframe2$ccmethod<-as.factor(dframe2$ccmethod)
dframe2$ccmethod = with(dframe2, factor(dframe2$ccmethod, levels = rev(levels(dframe2$ccmethod))))


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

dframe2$vi<-0
#res<-rma(yi=M,vi=M_sd^2,mods=~age+nsnps+samplesize+relatedness+ancestry,data=dframe2)
#res<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe+ancestry,data=dframe2)
#yi: vector of length k with the observed effect sizes or outcomes. See ‘Details’.
#vi: vector of length k with the corresponding sampling variances. See ‘Details’.
#res2<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+height+bmi+ccmethod+mean.probe+sd.probe,data=dframe2)
#res3<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe,data=dframe2)

w<-which(dframe2[,1]%in%c("MARS_omni"))
dframe3<-dframe2[-w,]
rma(yi=M,vi=vi,mods=~age,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt    0.0201  0.1063   0.1889  0.8501  -0.1882  0.2284   
#age       -0.0003  0.0024  -0.1393  0.8892  -0.0050  0.0043   


rma(yi=M,vi=vi,mods=~nmales,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt   -0.0089  0.1744  -0.0509  0.9594  -0.3506  0.3329   
#nmales     0.0347  0.3496   0.0994  0.9209  -0.6504  0.7199   


rma(yi=M,vi=vi,mods=~ancestry_ns,data=dframe3)
#Model Results:

#                  estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt             0.0216  0.0765   0.2827  0.7774  -0.1283  0.1715   
#ancestry_nssouth   -0.0307  0.1129  -0.2719  0.7857  -0.2519  0.1905   

rma(yi=M,vi=vi,mods=~nsnps,data=dframe3)

#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub    
#intrcpt    0.7567  0.2764   2.7377  0.0062   0.2150   1.2984  **
#nsnps     -0.0000  0.0000  -2.7579  0.0058  -0.0000  -0.0000  **

rma(yi=M,vi=vi,mods=~chip,data=dframe3)
#Model Results:

#          estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt     0.0221  0.0581   0.3803  0.7037  -0.0918  0.1359   
#chipEPIC   -0.1797  0.2040  -0.8809  0.3784  -0.5795  0.2201   

rma(yi=M,vi=vi,mods=~ncpg,data=dframe3)

#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt    0.1954  0.3035   0.6437  0.5197  -0.3995  0.7903   
#ncpg      -0.0000  0.0000  -0.6298  0.5288  -0.0000  0.0000   

rma(yi=M,vi=vi,mods=~samplesize,data=dframe3)

#Model Results:

#            estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt      -0.0709  0.0992  -0.7144  0.4750  -0.2654  0.1236   
#samplesize    0.0001  0.0001   0.9538  0.3402  -0.0001  0.0003   

rma(yi=M,vi=vi,mods=~relatedness,data=dframe3)
#Model Results:

#               estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt         -0.0085  0.0633  -0.1337  0.8936  -0.1326  0.1157   
#relatednessno    0.0740  0.1362   0.5430  0.5871  -0.1930  0.3409   

rma(yi=M,vi=vi,mods=~maf,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub   
#intrcpt   -1.9414  0.9289  -2.0901  0.0366  -3.7620  -0.1209  *
#maf       10.0554  4.7846   2.1016  0.0356   0.6778  19.4330  *

rma(yi=M,vi=vi,mods=~info,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt    1.1166  3.7583   0.2971  0.7664  -6.2496  8.4827   
#info      -1.1499  3.8963  -0.2951  0.7679  -8.7865  6.4867   

rma(yi=M,vi=vi,mods=~lambda4,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt    0.8702  1.7594   0.4946  0.6209  -2.5782  4.3186   
#lambda4   -0.8543  1.7449  -0.4896  0.6244  -4.2742  2.5657   

rma(yi=M,vi=vi,mods=~lambda16,data=dframe3)
#Model Results:

#          estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt     0.7645  0.4426   1.7274  0.0841  -0.1029  1.6319  .
#lambda16   -0.7778  0.4514  -1.7233  0.0848  -1.6625  0.1068  .

rma(yi=M,vi=vi,mods=~ccmethod,data=dframe3)
#Model Results:

#              estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt        -0.1947  0.1674  -1.1636  0.2446  -0.5227  0.1333   
#ccmethodcord    0.2268  0.1772   1.2797  0.2007  -0.1206  0.5741   

rma(yi=M,vi=vi,mods=~mean.probe,data=dframe3)
#Model Results:

#            estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt      -1.4451  1.2665  -1.1410  0.2539  -3.9275  1.0373   
#mean.probe    2.9440  2.5644   1.1480  0.2510  -2.0821  7.9700   

rma(yi=M,vi=vi,mods=~sd.probe,data=dframe3)
#Model Results:

#          estimate      se     zval    pval     ci.lb     ci.ub     
#intrcpt     1.0322  0.1779   5.8024  <.0001    0.6835    1.3808  ***
#sd.probe  -32.5863  5.5136  -5.9102  <.0001  -43.3927  -21.7798  ***

rma(yi=M,vi=vi,mods=~ancestry,data=dframe3)
#Model Results:

#                         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt                   -0.0196  0.0922  -0.2130  0.8313  -0.2004  0.1611   
#ancestryThe Netherlands    0.1788  0.1597   1.1193  0.2630  -0.1343  0.4918   
#ancestrySweden             0.4304  0.3325   1.2945  0.1955  -0.2212  1.0820   
#ancestrySpain             -0.1120  0.2062  -0.5433  0.5869  -0.5162  0.2921   
#ancestryScotland          -0.3338  0.2440  -1.3681  0.1713  -0.8119  0.1444   
#ancestryNew Zealand        0.4096  0.3325   1.2321  0.2179  -0.2420  1.0613   
#ancestryItaly             -0.0852  0.3325  -0.2562  0.7978  -0.7368  0.5665   
#ancestryGermany           -0.4867  0.3325  -1.4639  0.1432  -1.1383  0.1649   
#ancestryFrance             0.3192  0.3325   0.9602  0.3370  -0.3324  0.9709   
#ancestryFinland           -0.1482  0.2062  -0.7187  0.4723  -0.5523  0.2559   
#ancestryEstonia           -0.0312  0.2440  -0.1279  0.8983  -0.5094  0.4470   
#ancestryDenmark            0.0872  0.3325   0.2623  0.7931  -0.5644  0.7388   
#ancestryCanada            -0.2349  0.3325  -0.7064  0.4799  -0.8865  0.4168   
#ancestryAustralia          0.5016  0.2440   2.0562  0.0398   0.0235  0.9798  *

rma(yi=M,vi=vi,mods=~height,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt   -0.7836  0.4700  -1.6672  0.0955  -1.7048  0.1376  .
#height     0.4590  0.2866   1.6018  0.1092  -0.1026  1.0207   

rma(yi=M,vi=vi,mods=~bmi,data=dframe3)
#Model Results:

#         estimate      se     zval    pval    ci.lb    ci.ub    
#intrcpt    1.9751  0.6119   3.2277  0.0012   0.7757   3.1744  **
#bmi      -0.0800  0.0243  -3.2869  0.0010  -0.1277  -0.0323  **


#Model Results:

#                     estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt               -0.0398  0.1862  -0.2136  0.8309  -0.4046  0.3251   
#norm.method_fnFN      -0.0630  0.1996  -0.3154  0.7524  -0.4542  0.3283   
#norm.method_fnOther    0.2019  0.2065   0.9777  0.3282  -0.2029  0.6067   

rma(yi=M,vi=vi,mods=~norm.method,data=dframe3)
#Model Results:

#                      estimate      se     zval    pval    ci.lb   ci.ub   
#intrcpt                -0.0398  0.1945  -0.2044  0.8381  -0.4210  0.3415   
#norm.methodBMIQ         0.3426  0.3890   0.8806  0.3785  -0.4199  1.1051   
#norm.methodCPACOR      -0.0179  0.2751  -0.0652  0.9480  -0.5571  0.5212   
#norm.methodDasen        0.2515  0.2461   1.0221  0.3067  -0.2308  0.7337   
#norm.methodFN          -0.0630  0.2086  -0.3018  0.7628  -0.4718  0.3459   
#norm.methodglm          0.1547  0.3076   0.5030  0.6150  -0.4481  0.7575   
#norm.methodmethylumi    0.4298  0.3890   1.1047  0.2693  -0.3327  1.1923   
#norm.methodSWAN         0.3394  0.3890   0.8723  0.3830  -0.4231  1.1019   


rma(yi=M,vi=vi,mods=~nsnps+maf+bmi,data=dframe3)
rma(yi=M,vi=vi,mods=~nsnps+maf+bmi+ancestry,data=dframe3)



res<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe+ancestry,data=dframe3)
#yi: vector of length k with the observed effect sizes or outcomes. See ‘Details’.
#vi: vector of length k with the corresponding sampling variances. See ‘Details’.
res2<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+height+bmi+ccmethod+mean.probe+sd.probe,data=dframe3)
res3<-rma(yi=M,vi=vi,mods=~nsnps+chip+ncpg+samplesize+relatedness+maf+info+lambda16+lambda4+age+nmales+ccmethod+mean.probe+sd.probe,data=dframe3)

res2<-rma(yi=M,vi=vi,mods=~relatedness,data=dframe3)



p1<-ggplot(dframe3,aes(x=bmi,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average BMI") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="../images/MvaluebyBMI.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=maf,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average MAF") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="../images/Mvaluebymaf.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=nsnps,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("nSNPs") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="../images/Mvaluebynsnps.pdf",width=8,height=6)

p1<-ggplot(dframe3,aes(x=sd.probe,y=M))+
geom_point(aes(colour=factor(dframe3$study_names_in),size=dframe3$samplesize)) +
xlab("average SD across all probes") +
ylab("M Statistic") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="../images/Mvaluebysd.pdf",width=8,height=6)

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
pdf("mstat.forestplot.pdf",height=6,width=6)
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
