library(ggplot2)

load("../results/enrichments/snpcontrolsets_selection.rdata")
w<-which(is.na(f.all$snp_cis))
f.all$Category<-as.character(f.all$snp_cis)
f.all$Category[w]<-"no_mqtl"

f.all$Category<-gsub("TRUE","cisonly",f.all$Category)
f.all$Category<-gsub("FALSE","transonly",f.all$Category)
f.all$Category<-gsub("ambivalent","cis+trans",f.all$Category)

f.all$min_log10pval<-f.all$min_pval
w0<-which(f.all$min_pval==0)
mx<-min(f.all$min_pval[-w0],na.rm=T)
f.all$min_log10pval[w0]<-mx
f.all$min_log10pval<--log10(as.numeric(f.all$min_log10pval))

w<-which(f.all$mqtl_clumped=="TRUE")
f.all2<-f.all[w,]

r2<-read.table("./height_sds/vars.Height_Yengo_2018_mqtl_7cols.gz",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./height_sds/vars.Height_Yengo_2018_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,8]),":SNP")
f.all$height_mqtl<-"no mqtl"

w<-which(f.all$mqtl_clumped=="TRUE")
f.all$height_mqtl[w]<-"clumped mqtl"

w<-which(f.all$SNP%in%mqtl)
f.all$height_mqtl[w]<-"height mqtl"

sds<-paste0("chr",r3[,8],":SNP")
w<-which(f.all$SNP%in%sds)
f.all$height_mqtl[w]<-"height mqtl sds"
#ES = 2β^2f(1 − f)
f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

p1<-ggplot(f.all, aes(es, colour=height_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Mqtl_GeneticVariance.pdf")

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

#table(h_all$height_mqtl)

#    height mqtl height mqtl sds     height SNPs 
#             3948               44             653 

library(dplyr)
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

save(h_all,file="./height_sds/height_Yenko_plot.Robj")

