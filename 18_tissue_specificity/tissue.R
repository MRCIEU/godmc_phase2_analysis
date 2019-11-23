library(getmstatistic)  # for calculating M statistics
library(gridExtra)       # for generating tables
library(ggplot2)
library(data.table)
library(dplyr)

load("../results/16/16_clumped.rdata")
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
clumped$id<-paste0(clumped$snp,"_",clumped$cpg)
clumped$Allele1<-toupper(clumped$Allele1)
clumped$Allele2<-toupper(clumped$Allele2)
clumped$id2<-paste0(clumped$id,"_",clumped$Allele1,"_", clumped$Allele2)
data=as.data.table(clumped)
data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]
data[,cpg_cis:=ifelse(all(cis),"cis only",ifelse(all(!cis),"trans only","cis+trans")),by=c("cpgchr","cpgpos")]
clumped<-data
clumped$cpg_cis2<-clumped$cpg_cis
w<-which(clumped$cpg_cis=="cis+trans"&clumped$cis=="TRUE")
clumped$cpg_cis2[w]<-"cis+trans_cis"
w<-which(clumped$cpg_cis=="cis+trans"&clumped$cis=="FALSE")
clumped$cpg_cis2[w]<-"cis+trans_trans"


a<-read.table("adipose_16pairs.txt",he=T)
nrow(a)
table(a$EA)

m<-match(clumped$id,a$MARKERNAME)
a2<-data.frame(clumped,a[m,])

a2<-a2[which(!is.na(m)),]

a2$A1<-NA
a2$A2<-NA
i1 <- a2$EA==a2$Allele1&a2$NEA == a2$Allele2
i2 <- a2$NEA == a2$Allele1&a2$EA == a2$Allele2

a2$EA<-as.character(a2$EA)
a2$NEA<-as.character(a2$NEA)

a2$BETA2[i1]<-a2$BETA[i1]
a2$BETA2[i2]<--a2$BETA[i2]
a2$A1[i1]<-a2$EA[i1]
a2$A1[i2]<-a2$NEA[i2]
a2$A2[i1]<-a2$NEA[i1]
a2$A2[i2]<-a2$EA[i2]
a2$AF[i1]<-a2$EAF[i1]
a2$AF[i2]<-1-a2$EAF[i2]

table(i1,i2)
w<-which(i1=="FALSE"&i2=="FALSE")
if(length(w)>0){a2<-a2[-w,]}

p4<-ggplot(a2,aes(x=Effect,y=BETA2))+
geom_point() +
labs(x="GoDMC effect size",y="adipose effect size") +
facet_wrap(.~cpg_cis2,ncol=2,nrow=2)
ggsave(plot=p4, file="./images/adipose_blood.pdf", width=8, height=7)

gp = group_by(a2, cpg_cis2)
dplyr::summarize(gp, cor(Effect, BETA2))
# A tibble: 4 x 2
#  cpg_cis2        `cor(Effect, BETA2)`
#  <chr>                          <dbl>
#1 cis only                       0.700
#2 cis+trans_cis                  0.772
#3 cis+trans_trans                0.746
#4 trans only                     0.881

a3<-a2[which(a2$pval<1e-14),]
gp = group_by(a3, cpg_cis2)
dplyr::summarize(gp, cor(Effect, BETA2))
# A tibble: 4 x 2
#  cpg_cis2        `cor(Effect, BETA2)`
#  <chr>                          <dbl>
#1 cis only                       0.705
#2 cis+trans_cis                  0.774
#3 cis+trans_trans                0.746
#4 trans only                     0.881

a2<-data.frame(a2[,c("MARKERNAME","BETA2","SE","PVAL","N","AF")])

nrow(a)
nrow(a2)

w<-which(clumped$id%in%a2$id)
nomatch<-clumped[-w,"id"]
length(which(a$id%in%nomatch))

####
load("brainQTLs.rdat")
head(brainQTLs)
b<-brainQTLs

table(b$A1)
b$id<-paste0(b$SNP,"_",b$gene)
b$SE<-b$beta/b$t.stat
m<-match(clumped$id,b$id)
b2<-data.frame(clumped,b[m,])

b2<-b2[which(!is.na(m)),]
which(is.na(b2$beta))

b2$A1.orig<-b2$A1
b2$A2.orig<-b2$A2
b2$A1<-NA
b2$A2<-NA
i1 <- b2$A1.orig==b2$Allele1&b2$A2.orig == b2$Allele2
i2 <- b2$A2.orig == b2$Allele1&b2$A1.orig == b2$Allele2

b2$A1.orig<-as.character(b2$A1.orig)
b2$A2.orig<-as.character(b2$A2.orig)

b2$BETA2[i1]<-b2$beta[i1]
b2$BETA2[i2]<--b2$beta[i2]
b2$A1[i1]<-b2$A1.orig[i1]
b2$A1[i2]<-b2$A2.orig[i2]
b2$A2[i1]<-b2$A2.orig[i1]
b2$A2[i2]<-b2$A1.orig[i2]

table(i1,i2)
w<-which(i1=="FALSE"&i2=="FALSE")
if(length(w)>0){b2<-b2[-w,]}


p4<-ggplot(b2,aes(x=Effect,y=BETA2))+
geom_point() +
labs(x="GoDMC effect size",y="brain_effect size") +
facet_wrap(.~cpg_cis,ncol=2,nrow=2)
ggsave(plot=p4, file="./images/brain_blood.pdf", width=8, height=7)

gp = group_by(b2, cpg_cis2)
dplyr::summarize(gp, cor(Effect, BETA2))

#dplyr::summarize(gp, cor(Effect, BETA2))
# A tibble: 4 x 2
#  cpg_cis2        `cor(Effect, BETA2)`
#  <chr>                          <dbl>
#1 cis only                       0.456
#2 cis+trans_cis                  0.443
#3 cis+trans_trans                0.426
#4 trans only                     0.710

gp = group_by(b3, cpg_cis2)
dplyr::summarize(gp, cor(Effect, BETA2))
# A tibble: 4 x 2
#  cpg_cis2        `cor(Effect, BETA2)`
#  <chr>                          <dbl>
#1 cis only                       0.462
#2 cis+trans_cis                  0.448
#3 cis+trans_trans                0.426
#4 trans only                     0.710

b2<-data.frame(b2[,c("id","BETA2","SE","p.value")])

w<-which(clumped$id%in%b2$id)
nomatch<-clumped[-w,"id"]
length(which(b$id%in%nomatch))

nrow(b)
nrow(b2)
###

spl<-strsplit(b2$id,split="_")
spl<-do.call("rbind",spl)
b2<-data.frame(snp=spl[,1],cpg=spl[,2],b2)

h<-read.table("Hannon2016.txt",sep="\t",he=T)
w<-which(h$cpgchr==h$snpchr&abs(h$snppos-h$cpgpos)<1e6)
h2<-h[w,]


m<-match(h2$cpg,b2$cpg)
h2<-b2[m,]

save(clumped,a2,b2,h2,file="tissue.RData")

#
