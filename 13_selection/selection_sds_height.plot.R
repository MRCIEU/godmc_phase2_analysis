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

r2<-read.table("./height_sds/vars.GIANT_HEIGHT_mqtl",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./height_sds/vars.GIANT_HEIGHT_sds",sep=" ",he=T)
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

h<-read.table("GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
m<-match(h$MarkerName,bim$V2)
h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)
h$b_sq<-h$b^2
h$MAF<-h$Freq.Allele1.HapMapCEU
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
#r<-read.table("/panfs/panasas01/sscm/epzjlm/repo/goya/godmc/resources/genetics/wood.height.snps_af0.1.txt",he=F,stringsAsFactors=F)
r<-read.table("height_snps.txt",he=F,sep="\t")

h$height_mqtl<-"no_mqtl"
w<-which(h$MarkerName%in%r$V1)
h$height_mqtl[w]<-"height SNPs (n=697)"
h_all<-h[w,]

w<-which(h$SNP%in%mqtl)
h$height_mqtl[w]<-"height mqtl (n=39)"
h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$height_mqtl[w]<-"height mqtl sds (n=5)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)

#table(h_all$height_mqtl)

#    height mqtl height mqtl sds     height SNPs 
#             39               5             697 

p1<-ggplot(h_all, aes(es, colour=height_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Height_GeneticVariance.pdf")


