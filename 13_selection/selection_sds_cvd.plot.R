library(ggplot2)
library(dplyr)

#r2<-read.table("./cvd_sds/vars.Cardiogram_2_5cols",sep="\t",he=T)
r2<-read.table("./cvd_sds/vars.CAD_Nikpay_2015",sep=" ",he=T)
r2[,7]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./cvd_sds/cardiogram_sds_variants_6cols",sep="\t",he=T) #2
r3[,7]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,7]),":SNP")
sds<-paste0("chr",r3[,7],":SNP")

h<-read.table("./cvd_sds/cad.add.160614.website.txt",he=T)

table(as.character(h$effect_allele))

#      A       C       D       G       I       T 
#1927136 2380531  365339 2384400  470965 1927407 

h$INDEL<-"SNP"
w<-which(h$effect_allele%in%c("I","D"))
h$INDEL[w]<-"INDEL"
h$SNP<-paste0("chr",h$chr,":",h$bp_hg19,":",h$INDEL)

#bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
#m<-match(h$markername,bim$V2)
#h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)

h$b_sq<-h$beta^2
h$MAF<-h$effect_allele_freq
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))

r<-read.table("./cvd_sds/cvd_variants.txt",he=F,sep="\t")

h$cvd_mqtl<-"no_mqtl"
w<-which(h$markername%in%r$V1)
h$cvd_mqtl[w]<-"cvd SNPs (n=202)"
h_all<-h[w,]

w<-which(mqtl%in%h$SNP==F)
spl<-do.call("rbind",strsplit(mqtl[w],split=":"))
#[1] "chr2:204062727:SNP" "chr3:138099161:SNP" "chr6:160434842:SNP"
#[4] "chr2:203825427:SNP"

#h[h$bp_hg19%in%spl[,2],]
#              markername chr   bp_hg19 effect_allele noneffect_allele
#1382599 chr2:203825427:D   2 203825427             I                D
#1383040 chr2:204062727:D   2 204062727             I                D
#1988270 chr3:138099161:I   3 138099161             D                I
#4111028 chr6:160434842:I   6 160434842             D                I
#        effect_allele_freq median_info model      beta    se_dgc    p_dgc
#1382599           0.889278     0.99800 FIXED -0.132497 0.0155023 1.26e-17
#1383040           0.896366     0.97800 FIXED -0.135633 0.0162366 6.62e-17
#1988270           0.837203     0.99368 FIXED -0.074437 0.0125367 2.89e-09
#4111028           0.987163     0.91200 FIXED -0.335175 0.0507585 4.02e-11
#        het_pvalue n_studies INDEL                  SNP        b_sq      MAF
#1382599   0.422334        47 INDEL chr2:203825427:INDEL 0.017555455 0.110722
#1383040   0.395598        47 INDEL chr2:204062727:INDEL 0.018396311 0.103634
#1988270   0.800220        48 INDEL chr3:138099161:INDEL 0.005540867 0.162797
#4111028   0.369463        38 INDEL chr6:160434842:INDEL 0.112342281 0.012837
#                 es         cvd_mqtl
#1382599 0.003457113          no_mqtl
#1383040 0.003417814          no_mqtl
#1988270 0.001510375 cvd SNPs (n=202)
#4111028 0.002847250          no_mqtl


w<-which(h$SNP%in%mqtl)
h$cvd_mqtl[w]<-"cvd mqtl (n=15)"
#h$cvd_mqtl[w]<-"cvd mqtl (n=55)"
h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$cvd_mqtl[w]<-"cvd mqtl sds (n=2)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)

#table(h_all$cvd_mqtl)

#   cvd mqtl cvd mqtl sds     cvd SNPs 
#          15            2          202 

h_all%>%group_by(cvd_mqtl)%>%summarise(es=mean(es))
# A tibble: 3 x 2
#            cvd_mqtl          es
#               <chr>       <dbl>
#1    cvd mqtl (n=15) 0.001862760
#2 cvd mqtl sds (n=2) 0.010436141
#3   cvd SNPs (n=202) 0.001571966

#cvd_mqtl                es
#  <chr>                <dbl>
#1 cvd mqtl (n=56)    0.00256
#2 cvd mqtl sds (n=2) 0.0104 
#3 cvd SNPs (n=202)   0.00157

p1<-ggplot(h_all, aes(es, colour=cvd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="cvd_GeneticVariance.pdf")

h_all$cvd_mqtl <- factor(h_all$cvd_mqtl, levels = c("cvd SNPs (n=202)","cvd mqtl (n=15)", "cvd mqtl sds (n=2)"))
p1<-ggplot(h_all, aes(y=es, x=cvd_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="CVD_GeneticVariance_boxplot.pdf")

#no mapping
load("/newshared/godmc/database_files/snps.rdata")
h_all[which(h_all$SNP%in%out_df2$name==F),c("SNP","cvd_mqtl")]

length(mqtl) #59
h_all[which(h_all$SNP%in%mqtl==F),c("SNP","cvd_mqtl")]
length(unique(h_all[which(h_all$SNP%in%mqtl),c("SNP")])) #55
unique(h_all[which(h_all$SNP%in%mqtl),c("SNP")])

length(sds) #2
h_all[which(h_all$SNP%in%sds==F),c("SNP","cvd_mqtl")]
length(unique(h_all[which(h_all$SNP%in%sds),c("SNP")])) #2
unique(h_all[which(h_all$SNP%in%sds),c("SNP")])

save(h_all,file="./cvd_sds/cvd_plot.Robj")

