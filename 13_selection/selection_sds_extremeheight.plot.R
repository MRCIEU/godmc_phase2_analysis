library(ggplot2)
library(dplyr)

r2<-read.table("./height_sds/vars.EXTREME_HEIGHT_Berndt_2013",sep=" ",he=T)
#r2<-read.table("./height_sds/vars.EXTREME_HEIGHT_mqtl_6cols",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./height_sds/vars.EXTREME_HEIGHT_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))
mqtl<-paste0("chr",unique(r2[,8]),":SNP")

sds<-paste0("chr",r3[,8],":SNP")


#how much of the height variance

h<-read.table("./height_sds/GIANT_EXTREME_HEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt.gz",he=T)
load("/newshared/godmc/database_files/snps.rdata")


table(h$Allele1)

#     a      c      g      t 
#482299 502872 501689 479696 

#bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
#m<-match(h$MarkerName,bim$V2)

m<-match(h$MarkerName,out_df2$rsid)
h<-data.frame(SNP=out_df2[m,"name"],h)
length(which(is.na(h$SNP))) #14780

h$b_sq<-h$b^2
h$MAF<-h$Allele1_Freq_HapMapCEU
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
#r<-read.table("/panfs/panasas01/sscm/epzjlm/repo/goya/godmc/resources/genetics/wood.height.snps_af0.1.txt",he=F,stringsAsFactors=F)
load("/newshared/godmc/database_files/snps.rdata")

m<-match(h$MarkerName,out_df2$rsid)
id<-out_df2[m,c("name","rsid")]
h<-data.frame(id,h)

#r<-read.table("./height_sds/extreme_height_snps.txt",he=F,sep="\t")
r<-read.table("./height_sds/ext_height_rsid",he=F)
h$height_mqtl<-"no_mqtl"
#w<-which(h$MarkerName%in%r$V1)
w<-which(h$MarkerName%in%r$V1)
h$height_mqtl[w]<-rep("extreme height SNPs (n=60)",length(w))
h_all<-h[w,]

w<-which(h$name%in%mqtl)
h$height_mqtl[w]<-"extreme height mqtl (n=48)"
h_all2<-h[w,]

w<-which(h$name%in%sds)
h$height_mqtl[w]<-"extreme height mqtl sds (n=4)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)


table(h_all$height_mqtl)

#extreme height mqtl sds (n=4)    extreme height SNPs (n=60) 
#                            8                            60 
#                  height mqtl 
#                           38 

p1<-ggplot(h_all, aes(es, colour=height_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="ExtremeHeight_GeneticVariance.pdf")

h_all$height_mqtl <- factor(h_all$height_mqtl, levels = c("extreme height SNPs (n=60)","extreme height mqtl (n=48)", "extreme height mqtl sds (n=4)"))
p1<-ggplot(h_all, aes(y=es, x=height_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="Extremeheight_GeneticVariance_boxplot.pdf")



h_all%>%group_by(height_mqtl)%>%summarise(es=mean(es,na.rm=T))

# A tibble: 3 x 2
#                    height_mqtl          es
#                          <chr>       <dbl>
#1    extreme height mqtl (n=48) 0.015082114
#2 extreme height mqtl sds (n=4) 0.013163220
#3    extreme height SNPs (n=60) 0.005124636

save(h_all,file="./height_sds/extremeheight_plot.Robj")

