library(getmstatistic)  # for calculating M statistics
library(gridExtra)       # for generating tables
library(ggplot2)
source("plot.R")

#getmstatistic uses summary data i.e. study effect-sizes and their corresponding standard errors
#to calculate M statistics (One M for each study in the meta-analysis).
#In particular, getmstatistic employs the inverse-variance weighted random effects regression
#model provided in the metafor R package to extract SPREs (standardized predicted random effects)
#which are then aggregated to formulate M statistics.

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
#[1] 337
plots_dir <-"./images/pdf"
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
#getmstatistic_results

#filename_mstats_vs_avg_effectsize <- base::paste0("Mstatistics_vs_average_variant_effectsize_", nstudies, "studies_", nsnps, "snps.tif")
#grDevices::tiff(filename_mstats_vs_avg_effectsize, width = 23.35, height = 17.35, units = "cm", res = 300, compression = "lzw", pointsize = 14)
#h <- ggplot2::ggplot(usta_mean_scatter_strength, ggplot2::aes(base::log(oddsratio), usta_mean, colour = usta_mean, label = study_names)) + ggplot2::geom_point(size = 4.5) + ggplot2::geom_text(ggplot2::aes(label = study), hjust = 1.2, vjust = -0.5, size = 2.5, colour = "azure4") + ggplot2::scale_colour_gradientn(name = "M statistic", colours = grDevices::rainbow(11)) + ggplot2::scale_x_continuous(trans="log", limits=c(x_axis_min, x_axis_max), breaks=base::round(base::seq(x_axis_min, x_axis_max, x_axis_increment_in),x_axis_round_in), minor_breaks=ggplot2::waiver(), labels = base::round(base::exp(base::seq(x_axis_min, x_axis_max, x_axis_increment_in)),x_axis_round_in)) + ggplot2::theme_bw() + ggplot2::scale_fill_hue(c = 45, l = 40) + ggplot2::xlab("Average effect size (oddsratio)") + ggplot2::ylab("M statistic") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 14)) + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14), axis.text.y = ggplot2::element_text(size = 14))
#hi <- h + ggplot2::geom_hline(ggplot2::aes(yintercept = yval), data = dat_hlines, colour = "grey80", linetype = "dashed", lwd = 0.4) + ggplot2::theme(legend.text = ggplot2::element_text(size = 10)) + ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) + ggplot2::theme(legend.position = "bottom")
#base::print(hi + ggplot2::geom_hline(ggplot2::aes(yintercept = c(0,0)), data = dat_hlines, colour = "grey80", linetype = "solid", lwd = 0.4))
#grDevices::dev.off()

dframe$study_names_in<-gsub("00_ARIES","ARIES",dframe$study_names_in)
dframe$study_names_in<-gsub("as_cc","EGC_asthma",dframe$study_names_in)
dframe$study_names_in<-gsub("ccg","EGC_CTG",dframe$study_names_in)
dframe$study_names_in<-gsub("InterAct","EPIC_Norfolk",dframe$study_names_in)
dframe$study_names_in<-gsub("Leiden_Longevity_Study","LLS",dframe$study_names_in)
dframe$study_names_in<-gsub("Project_MinE_s27","MinE",dframe$study_names_in)
dframe$study_names_in<-gsub("IOW3g","IOW F2",dframe$study_names_in)

dframe$study_names_in<-paste0(dframe$study_names_in," (n=",dframe$beta_n,")")
p1 <- ggplot(dframe, aes(x=as.factor(study_names_in), y=beta_in)) +
geom_boxplot() +
theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_text(size=10)) +
labs(x="Study", y="Effect size")
ggsave(plot=p1, file="./images/effectsizechr20bystudy.v2.pdf", width=7, height=7)

library(plyr)
tmp<-daply(df, .(variant_names_in,study_names_in), function(x) x$beta_in)
tmp<-data.frame(tmp)

p<-ggpairs(tmp[,1:9])
ggsave(p,file="./images/betacomptest1.png")

p<-ggpairs(tmp[,10:18])
ggsave(p,file="./images/betacomptest2.png")

p<-ggpairs(tmp[,19:27])
ggsave(p,file="./images/betacomptest3.png")

p<-ggpairs(tmp[,28:36])
ggsave(p,file="./images/betacomptest4.png")

p<-ggpairs(tmp)
ggsave(p,file="./images/betacomptest.png",height=12,width=12)


#p<-ggpairs(
#  tmp,
 # upper = list(continuous = "cor"),
  #lower = list(continuous = wrap("points",size=0.2)),
  #diag = list(continuous = wrap('diagAxis', labelSize = 2,gridLabelSize=0)),
  #columnLabels=rep("",ncol(tmp)),
#)



