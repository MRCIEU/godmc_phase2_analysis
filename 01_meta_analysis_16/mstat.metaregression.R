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

# Sort getmstatistic_results dataframe by M statistics



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






