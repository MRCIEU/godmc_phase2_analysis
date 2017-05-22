
#extract results files
#results_dir="/projects/MRC-IEU/groups/godmc/meta-analysis/input"
#metal_in="16_${i}.in"
#cd $results_dir

#for cohort in ${cohort_dir}*16.tar
#do

#	echo $cohort
#	cohortname=$(basename "$cohort" .tar)
#	echo $cohortname
#	mkdir -p $results_dir/${cohortname}
#	cd $results_dir/${cohortname}
#	tar -xvf ../${cohortname}.tar results/16/results_*.gz
#	cd $results_dir
	
#done




#cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16

cohort_dir="/projects/MRC-IEU/groups/godmc/meta-analysis/input"
result_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/mstat"

metal_dir="../scratch/16_${i}"
cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
library(getmstatistic)  # for calculating M statistics
library(gridExtra)       # for generating tables

head(heartgenes214)
#  beta_flipped      gcse          variants studies cases controls fdr214_gwas46
#1  -0.06325500 0.1607954  chr10:75595440:D      01   278      312             1
#2  -0.00041045 0.1443252 chr11:100634736:D      01   278      312             1
#3   0.07529900 0.1415288  chr11:77185247:I      01   278      312             1
#4   0.06133000 0.1893356  chr12:76386403:I      01   278      312             1
#5   0.19341000 0.1362752  chr14:75614504:I      01   278      312             1
#6  -0.03793500 0.1312889  chr16:75308440:D      01   278      312             1

path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis"

load(paste0(path,"/results/16/16_clumped.rdata"))
dim(clumped)
#[1] 312024     21
a<-clumped
a$id<-as.character(paste(a$snp,a$cpg,sep="_"))
a14.cis.out<-a[which(a$cis==F & a$pval<1e-14),]
dim(a14.cis.out)
#[1] 156857     22
o<-order(a14.cis.out$pval)
a14.cis.out<-a14.cis.out[o,]
probe<-unique(a14.cis.out$cpg)
m<-match(probe,a14.cis.out$cpg)
a14.cis.out<-a14.cis.out[m,]
#a14.cis.out<-a[which(a$cis==T & a$pval<1e-14&a$cpgchr=="chr20"),]
#dim(a14.cis.out)
#[1] 3213   21
#dim(data.frame(table(a14.cis.out$snp)))
#[1] 2497    2

o<-order(as.numeric(as.character(a14.cis.out$cpgpos)))
a14.cis.out<-a14.cis.out[o,]

l<-list.files(path=cohort_dir)
l2<-list.files(paste(cohort_dir,"/",l[1],"/results/16/",sep=""))
study<-gsub("_16","",l)

res<-data.frame()
for (i in 1:length(l)){
#for (i in 1:2){
cat(i,"\n")

	for (j in (1:25)){
		#for (j in 1:length(l2)){
		cat(j,"\n")
r<-read.table(paste(cohort_dir,"/",l[i],"/results/16/results_",j,".gz",sep=""),he=T)
m<-which(r$SNP%in%a14.cis.out$snp)

if(length(m)>0){

r2<-data.frame(r[m,],study=study[i])
res<-rbind(res,r2)
}
}
}

res$study<-as.character(res$study)
res$SNP<-as.character(res$SNP)
res$MARKERNAME<-as.character(res$MARKERNAME)

res2<-res[which(res$PVAL<1e-14),]

m<-data.frame(table(res2$MARKERNAME))
m<-m[which(m$Freq>1),1]
res2<-res2[res2$MARKERNAME%in%m,]

length(unique(res2$MARKERNAME))
#[1] 22813
getmstatistic_results <- getmstatistic(res2$BETA, res2$SE,res2$MARKERNAME, res2$study)

dframe <- getmstatistic_results$M_dataset
head(dframe)

# Retrieve dataset of stronger than average studies (significant at 5% level)
getmstatistic_results$influential_studies_0_05
 
# Retrieve dataset of weaker than average studies (significant at 5% level)
getmstatistic_results$weaker_studies_0_05
 
# Retrieve number of studies and variants
getmstatistic_results$number_studies
getmstatistic_results$number_variants
 
# Retrieve expected mean, sd and critical M value at 5% significance level
getmstatistic_results$M_expected_mean
getmstatistic_results$M_expected_sd
getmstatistic_results$M_crit_alpha_0_05

save.image("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/mstat/mstats.RData")
