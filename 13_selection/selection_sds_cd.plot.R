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

#r2<-read.table("./cd_sds/vars.GIANT_HEIGHT_mqtl",sep=" ",he=T)
#r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./cd_sds/vars.IBD_CD_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

#mqtl<-paste0("chr",unique(r2[,8]),":SNP")

f.all$cd_mqtl<-"no mqtl"

w<-which(f.all$mqtl_clumped=="TRUE")
f.all$cd_mqtl[w]<-"clumped mqtl"

#w<-which(f.all$SNP%in%mqtl)
#f.all$cd_mqtl[w]<-"cd mqtl"

sds<-paste0("chr",r3[,8],":SNP")
w<-which(f.all$SNP%in%sds)
w<-which(sds%in%f.all$SNP==F)
if(length(w)>0){
sds[w]<-gsub(":SNP",":INDEL",sds[w])}

w<-which(f.all$SNP%in%sds)
f.all$cd_mqtl[w]<-"cd mqtl sds"

#ES = 2β^2f(1 − f)
f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

p1<-ggplot(f.all, aes(es, colour=cd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Mqtl_GeneticVariance.pdf")

#how much of the cd variance

h<-read.table("./cd_sds/EUR.CD.gwas_info03_filtered.assoc",he=T)
w1<-which(h$A1=="D")
w2<-which(h$A2=="D")
w<-unique(c(w1,w2))
h$SNP<-as.character(h$SNP)
h$SNP<-paste0("chr",h$CHR,":",h$BP,":","SNP")
h[w,"SNP"]<-paste0("chr",h[w,"CHR"],":",h$BP[w],":","INDEL")
#bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
#m<-match(h$SNP,bim$V2)
#h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)
h$b_sq<-h$OR^2
h$MAF<-h$FRQ_A_5956
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))

r<-read.table("./cd_sds/cd_snps.txt",he=F,sep="\t")

h$cd_mqtl<-"no_mqtl"
w<-which(h$SNP.1%in%r$V1)
h$cd_mqtl[w]<-"cd SNPs (n=33)"
h_all<-h[w,]

#w<-which(h$SNP%in%mqtl)
#h$cd_mqtl[w]<-"cd mqtl"
#h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$cd_mqtl[w]<-"cd mqtl sds (n=3)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)
h_all<-rbind(h_all,h_all3)

table(h_all$cd_mqtl)

#cd mqtl sds (n=3)    cd SNPs (n=33) 
#                3                33 

p1<-ggplot(h_all, aes(es, colour=cd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="cd_GeneticVariance.pdf")

library(dplyr)
h_all%>%group_by(cd_mqtl)%>%summarise(es=mean(es,na.rm=T))

# A tibble: 2 x 2
#            height_mqtl          es
#                  <chr>       <dbl>
#1 height mqtl sds (n=4) 0.013163220
#2    height SNPs (n=60) 0.005124636



