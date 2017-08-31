#extract results files
#/panfs/panasas01/shared-godmc/sftp/GoDMC/ecarnero/TwinsUK/results/01
cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"


cd $results_dir

for cohort in ${cohort_dir}*01.tgz
do
echo $cohort


	cohortname=$(basename "$cohort" .tgz)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}
	cd $cohort_dir/${cohortname}


	#tar -xvf ../${cohortname}.tgz results/01/cohort_descriptives.RData
    #tar -xvf ../${cohortname}.tgz config
    #related=`grep "related=" config |perl -pe 's/related=//g'|perl -pe 's/"//g'`
    #echo $cohortname $related >>/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt 
    tar -xvf ../${cohortname}.tgz results/01/methylation_summary.RData
done
###
library(ggplot2)
library(metafor)
library(dplyr)   # for sorting data.frames
library(getmstatistic)  # for calculating M statistics
library(gridExtra)

cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/"

load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/mstat/mstats.RData")
dframe2<-unique(data.frame(study_names_in=dframe$study_names_in,M=dframe$M,M_sd=dframe$M_sd,M_se=dframe$M_se))
dframe2[,1]<-gsub("00_ARIES","ARIES",dframe2[,1])
path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis"
w<-which(dframe2[,1]%in%c("DunedinAge38"))
dframe2<-dframe2[-w,]



dframe2$age<-NA
dframe2$nsnps<-NA
dframe2$relatedness<-NA
for (i in 1:nrow(dframe2)){
cat(dframe2[i,1],"\n")
load(paste0(cohort_dir,dframe2[i,1],"_01/results/01/cohort_descriptives.RData"))
age<-cohort_summary$mqtl_mean_age
dframe2$age[i]<-age
nsnps<-cohort_summary$n_snp
dframe2$nsnps[i]<-nsnps
#con<-read.table(paste0(cohort_dir,dframe2[i,1],"_01/config"),he=F,sep=" ",stringsAsFactors=F,skip=1)
#con<-scan(paste0(cohort_dir,dframe2[i,1],"_01/config"),skip=1)
#g<-grep("related=",con$V1)
#rel<-con$V1[g]
#rel<-gsub("related=","",rel)
#dframe2$relatedness[i]<-rel
}

ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizes.txt")
names(ss)<-c("study","n")
ss[,1]<-gsub("_16","",ss[,1])
ss[,1]<-gsub("00_ARIES","ARIES",ss[,1])
m<-match(dframe2$study_names_in,ss[,1])
dframe2$samplesize<-ss[m,-1]

rel<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt")
rel[,1]<-gsub("_01","",rel[,1])
m<-match(dframe2$study_names_in,rel[,1])
dframe2$relatedness<-rel[m,-1]

ancestry<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortancestry.txt",sep="\t")
m<-match(dframe2$study_names_in,ancestry[,1])
dframe2$ancestry<-ancestry[m,-1]

dframe2$vi<-0
#res<-rma(yi=M,vi=M_sd^2,mods=~age+nsnps+samplesize+relatedness+ancestry,data=dframe2)
res<-rma(yi=M,vi=vi,mods=~age+nsnps+samplesize+relatedness+ancestry,data=dframe2)
#yi: vector of length k with the observed effect sizes or outcomes. See ‘Details’.
#vi: vector of length k with the corresponding sampling variances. See ‘Details’.

w<-which(dframe2[,1]%in%c("PRECISESADS"))
dframe3<-dframe2[-w,]
res<-rma(yi=M,vi=vi,mods=~age+nsnps+samplesize+relatedness+ancestry,data=dframe3)

p1<-ggplot(dframe2,aes(x=samplesize,y=nsnps))+
geom_point(aes(colour=factor(dframe2$study_names_in),size=dframe2$samplesize)) +
xlab("Cohort N") +
ylab("N SNPs") +
theme(axis.text.x = element_text(face = "bold"))
ggsave(p1,file="../images/SNPsbyNcohort.pdf",width=8,height=6)


# Sort getmstatistic_results dataframe by M statistics
getm_res_srtd <- dplyr::arrange(dframe, M)

# First, drop duplicate study_names in the sorted getmsatistic_results
getm_res_srtd_nodups <- subset(getm_res_srtd, getm_res_srtd$snp == 5)

# Checking dimensions to confirm that we have 48 studies
str(getm_res_srtd_nodups)


getm_res_plus_n <- merge(getm_res_srtd_nodups, ss, by.x="study_names_in", by.y="study")

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
text(1.6, 27.3, "M statistic [95% CI]", pos=2, cex=0.66)
text(-1.4, 27.3, "N", pos=4, cex=0.66)
#text(-1.39, 50, "Controls", pos=4, cex=0.66)
text(-3, 27.3, "Study", pos=4, cex=0.66)
dev.off()
