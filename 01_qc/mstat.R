library(getmstatistic)  # for calculating M statistics
library(gridExtra)       # for generating tables
library(ggplot2)

cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
result_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/mstat"

l<-list.files(path=cohort_dir)
w<-which(l%in%c("ARIES_16","MARS_omni_16","Factor_V_Leiden_Family_Study_16"))

if(length(w)>0){
l<-l[-w]}
study<-gsub("_16","",l)

cat(length(l),"cohorts","\n")
chunks<-read.table("chr20_chunks.txt",sep="\t")

res<-data.frame()
for (i in 1:length(l)){
cat(i,"\n")

for (chunk in 1:length(chunks[,1])){
j<-chunks[chunk,1]
cat(j,"\n")
p<-paste0("./chr20.assoc/",l[i],".chr20.assoc_",j,".txt")
f<-file.size(p)

if(f>0){
r<-read.table(p,he=F)
r<-data.frame(r,study=study[i])
res<-rbind(res,r)
}

}
}
names(res)<-c("MARKERNAME","PHENOTYPE","CISTRANS","BETA","SE","PVAL","N","EA","NEA","EAF","SNP","study")

cat(length(unique(res$MARKERNAME)),"\n")

res$study<-as.character(res$study)
res$SNP<-as.character(res$SNP)
res$MARKERNAME<-as.character(res$MARKERNAME)

#res2<-res[which(res$PVAL<1e-14),]

res2<-res
m<-data.frame(table(res2$MARKERNAME))
m<-m[which(m$Freq>1),1]
res2<-res2[res2$MARKERNAME%in%m,]

length(unique(res2$MARKERNAME))
#[1] 332
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

save(getmstatistic_results,dframe,file="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/mstat/mstats_chr20.RData")

p1 <- ggplot(dframe, aes(x=as.factor(study_names_in), y=beta_in)) +
geom_boxplot() +
theme(axis.title.x=element_blank(),axis.title.y=element_text(size=8),axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_text(size=6)) +
labs(x="study", y="effect size")
ggsave(plot=p1, file="./images/effectsizechr20bystudy.pdf", width=7, height=7)


