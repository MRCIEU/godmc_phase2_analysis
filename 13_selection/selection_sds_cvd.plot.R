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

r2<-read.table("./cvd_sds/vars.Cardiogram_2_5cols",sep="\t",he=T)
r2[,7]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./cvd_sds/cardiogram_sds_variants_6cols",sep="\t",he=T) #2
r3[,7]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,7]),":SNP")
f.all$cvd_mqtl<-"no mqtl"

w<-which(f.all$mqtl_clumped=="TRUE")
f.all$cvd_mqtl[w]<-"clumped mqtl"

w<-which(f.all$SNP%in%mqtl)
f.all$cvd_mqtl[w]<-"cvd mqtl (n=15)"

sds<-paste0("chr",r3[,7],":SNP")
w<-which(f.all$SNP%in%sds)
f.all$cvd_mqtl[w]<-"cvd mqtl sds (n=2)"
#ES = 2β^2f(1 − f)
f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

p1<-ggplot(f.all, aes(es, colour=cvd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Mqtl_GeneticVariance.pdf")

#how much of the height variance

h<-read.table("./cvd_sds/cad.add.160614.website.txt",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
bim<-which(bim$V2%in%h$markername)
m<-match(h$markername,bim$V2)
h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)
h$b_sq<-h$beta^2
h$MAF<-h$effect_allele_freq
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
#r<-read.table("/panfs/panasas01/sscm/epzjlm/repo/goya/godmc/resources/genetics/wood.height.snps_af0.1.txt",he=F,stringsAsFactors=F)
r<-read.table("./cvd_sds/cvd_variants.txt",he=F,sep="\t")

h$cvd_mqtl<-"no_mqtl"
w<-which(h$markername%in%r$V1)
h$cvd_mqtl[w]<-"cvd SNPs (n=202)"
h_all<-h[w,]

w<-which(h$SNP%in%mqtl)
h$cvd_mqtl[w]<-"cvd mqtl (n=15)"
h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$cvd_mqtl[w]<-"cvd mqtl sds (n=2)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)

#table(h_all$cvd_mqtl)

#   cvd mqtl cvd mqtl sds     cvd SNPs 
#          15            2          202 
library(dplyr)
h_all%>%group_by(cvd_mqtl)%>%summarise(es=mean(es))
# A tibble: 3 x 2
#            cvd_mqtl          es
#               <chr>       <dbl>
#1    cvd mqtl (n=15) 0.001862760
#2 cvd mqtl sds (n=2) 0.010436141
#3   cvd SNPs (n=202) 0.001571966


p1<-ggplot(h_all, aes(es, colour=cvd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="cvd_GeneticVariance.pdf")


