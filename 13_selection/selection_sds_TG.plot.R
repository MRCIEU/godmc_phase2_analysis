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

r2<-read.table("./TG_sds/vars.TG_6cols",sep=" ",he=T)
r2[,7]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #277

r3<-read.table("./TG_sds/vars.TG_sds",sep=" ",he=T) #2
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,7]),":SNP")
f.all$TG_mqtl<-"no mqtl"

w<-which(f.all$mqtl_clumped=="TRUE")
f.all$TG_mqtl[w]<-"clumped mqtl"

w<-which(f.all$SNP%in%mqtl)
f.all$TG_mqtl[w]<-"TG mqtl (n=117)"

sds<-paste0("chr",r3[,8],":SNP")
w<-which(f.all$SNP%in%sds)
f.all$TG_mqtl[w]<-"TG mqtl sds (n=4)"
#ES = 2β^2f(1 − f)
f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

p1<-ggplot(f.all, aes(es, colour=TG_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Mqtl_GeneticVariance.pdf")

#how much of the height variance

h<-read.table("./TG_sds/jointGwasMc_TG.txt.gz",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
bim<-data.frame(SNP=paste0("chr",bim[,1],":",bim[,4],":SNP"),bim)
h<-data.frame(SNP=paste0(h$SNP_hg19,":SNP"),h)
w<-which(bim$SNP%in%h$SNP)
bim<-bim[w,]
m<-match(h$SNP,bim$SNP)

#bim<-which(bim$V2%in%h$markername)
#m<-match(h$markername,bim$V2)
#h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)


h$b_sq<-h$beta^2
h$MAF<-h$Freq.A1.1000G.EUR
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))

r<-read.table("./TG_sds/Willer_Teslovich.txt",he=T,sep="\t")
r<-r[which(r$Trait=="TG"),]
#r<-r[r$Source=="Willer",]

h$TG_mqtl<-"no_mqtl"
w<-which(h$rsid%in%unique(r$SNPTEST_ID))
h$TG_mqtl[w]<-"TG SNPs (n=44)"
h_all<-h[w,]

w<-which(h$SNP%in%unique(mqtl))
h$TG_mqtl[w]<-"TG mqtl (n=117)"
h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$TG_mqtl[w]<-"TG mqtl sds (n=4)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)

#table(h_all$TG_mqtl)

#TG mqtl (n=116) TG mqtl sds (n=3)    TG SNPs (n=44) 
#              117                 4                44 

p1<-ggplot(h_all, aes(es, colour=TG_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="TG_GeneticVariance.pdf")

library(dplyr)
h_all%>%group_by(TG_mqtl)%>%summarise(es=mean(es,na.rm=T))
# A tibble: 3 x 2
#            TG_mqtl           es
#              <chr>        <dbl>
#1   TG mqtl (n=116) 0.0007856486
#2 TG mqtl sds (n=4) 0.0003527631
#3    TG SNPs (n=44) 0.0008762782









