library(ggplot2)
library(dplyr)

r2<-read.table("./height_sds/vars.Height_Yengo_2018",sep=" ",he=T)
#r2<-read.table("./height_sds/vars.Height_Yengo_2018_mqtl_7cols.gz",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./height_sds/vars.Height_Yengo_2018_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,8]),":SNP")

sds<-paste0("chr",r3[,8],":SNP")
#ES = 2β^2f(1 − f)

#how much of the height variance

h<-read.table("./height_sds/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz",he=T)
h<-data.frame(ID=paste0("chr",h$CHR,":",h$POS,":","SNP"),h)
a1<-(nchar(as.character(h$Tested_Allele)))
a2<-(nchar(as.character(h$Other_Allele)))
w1<-which(a1>1)
w2<-which(a2>1)
length(unique(c(w1,w2)))
w<-unique(c(w1,w2))
id<-paste0("chr",h$CHR,":",h$POS,":","INDEL")
h$ID<-as.character(h$ID)
h$ID[w]<-id

h$b_sq<-h$BETA^2
h$MAF<-h$Freq_Tested_Allele_in_HRS
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
#r<-read.table("/panfs/panasas01/sscm/epzjlm/repo/goya/godmc/resources/genetics/wood.height.snps_af0.1.txt",he=F,stringsAsFactors=F)

r<-read.table("./height_sds/Meta-analysis_Wood_et_al+UKBiobank_2018_top_3290_from_COJO_analysis.txt.gz",he=T)
h$height_mqtl<-"no_mqtl"
w<-which(h$SNP%in%r$SNP)
h$height_mqtl[w]<-"height SNPs (n=3290)"
h_all<-h[w,]

w<-which(h$ID%in%mqtl)
h$height_mqtl[w]<-"height mQTL(n=3948)"
h_all2<-h[w,]

w<-which(h$ID%in%sds)
h$height_mqtl[w]<-"height mQTL SDS (n=44)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)

table(h_all$height_mqtl)

#    height mqtl height mqtl sds     height SNPs 
#             3948               44             653 

h_all%>%group_by(height_mqtl)%>%summarise(es=mean(es))

p1<-ggplot(h_all, aes(es, colour=height_mqtl)) +
geom_density() +
xlim(0,0.001) +
labs(x="Genetic Variance")
ggsave(p1,file="Height_Yenko_GeneticVariance.pdf")

h_all$height_mqtl <- factor(h_all$height_mqtl, levels = c("height SNPs (n=3290)","height mQTL(n=3948)", "height mQTL SDS (n=44)"))
p1<-ggplot(h_all, aes(y=es, x=height_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="Height_Yenko_GeneticVariance_boxplot.pdf")

load("/newshared/godmc/database_files/snps.rdata")
h_all[which(h_all$SNP%in%out_df2$name==F),c("SNP","cd_mqtl")]
#                      SNP        cd_mqtl
#192126  chr10:35542343:SNP cd mqtl (n=60)
#3708494 chr16:50756540:SNP cd mqtl (n=60)
#8699915  chr6:32767249:SNP cd mqtl (n=60)

length(mqtl) #5033
h_all[which(h_all$ID%in%mqtl==F),c("SNP","height_mqtl")]
length(unique(h_all[which(h_all$ID%in%mqtl),c("SNP")])) #5033
unique(h_all[which(h_all$SNP%in%mqtl),c("SNP")])

length(sds) #44
h_all[which(h_all$ID%in%sds==F),c("SNP","height_mqtl")]
length(unique(h_all[which(h_all$ID%in%sds),c("SNP")])) #44
unique(h_all[which(h_all$ID%in%sds),c("SNP")])

save(h_all,file="./height_sds/height_Yenko_plot.Robj")

