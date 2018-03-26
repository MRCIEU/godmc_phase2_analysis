#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)

r<-read.table(paste0(path,"mqtl_tfbs/garfield.test.mqtl_tfbs.out"),he=T)
r$cis_snp<-"All"

cis<-read.table(paste0(path,"mqtl_cis_tfbs/garfield.test.mqtl_cis_tfbs.out"),he=T)
cis$cis_snp<-"cis only"

trans<-read.table(paste0(path,"mqtl_trans_tfbs/garfield.test.mqtl_trans_tfbs.out"),he=T)
trans$cis_snp<-"trans only"

amb<-read.table(paste0(path,"mqtl_ambivalent_tfbs/garfield.test.mqtl_ambivalent_tfbs.out"),he=T)
amb$cis_snp<-"ambivalent"

padj<-read.table(paste0(path,"mqtl_tfbs/garfield.Meff.mqtl_tfbs.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

df<-rbind(r,cis,trans,amb)
df$logOddsRatio<-log(df$OR)
df2<-df[df$PThresh=="1e-14",]
nvar<-df2%>%group_by(cis_snp)%>%summarize(mean(NThresh))

nvar
# A tibble: 4 x 2
#     cis_snp `mean(NThresh)`
#       <chr>           <dbl>
#1        All           55747
#2 ambivalent           10954
#3   cis only           52700
#4 trans only             928


df2$cis_snp<-gsub("All",paste0("All (N=",nvar[1,2]," SNPs)"),df2$cis_snp)
df2$cis_snp<-gsub("ambivalent",paste0("ambivalent (N=",nvar[2,2]," SNPs)"),df2$cis_snp)
df2$cis_snp<-gsub("cis only",paste0("cis only (N=",nvar[3,2]," SNPs)"),df2$cis_snp)
df2$cis_snp<-gsub("trans only",paste0("trans only (N=",nvar[4,2]," SNPs)"),df2$cis_snp)

df2$Type<-as.factor(df2$Type)
df2$Category<-as.factor(df2$Category)
cats<-levels(df2$Category)

pval_lim<-padj

w<-which(df2$Pvalue==0)
m<-min(df2$Pvalue[-w])
df2$Pvalue[w]<-m
df2$Type<-sub("_[^_]+$", "", df2$Type)

for (i in 1:(length(cats))){
cat(i,"\n")

if(cats[i]%in%c("TFBS")) {
df3<-df2[which(df2$Category==cats[i]&df2$cis_snp!=paste0("All (N=",nvar[1,2]," SNPs)")),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
facet_wrap(~cis_snp,scale="free_y",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/epigenetic_",cats[i],".pdf"),height=10,width=16)
}
}



group_by(df3[which(df3$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),]) %>% summarize(m = max(Pvalue))
group_by(df3[which(df3$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),]) %>% summarize(m = max(Pvalue))
group_by(df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),]) %>% summarize(m = max(Pvalue))

group_by(df3[which(df3$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),]) %>% summarize(m = min(Pvalue))
group_by(df3[which(df3$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),]) %>% summarize(m = min(Pvalue))
group_by(df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),]) %>% summarize(m = min(Pvalue))

group_by(df3[which(df3$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),],Tissue) %>% summarize(m = max(OR))
group_by(df3[which(df3$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),], Tissue) %>% summarize(m = max(OR))
group_by(df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),], Tissue) %>% summarize(m = max(OR))

group_by(df3[which(df3$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),],Tissue) %>% summarize(m = min(OR))
group_by(df3[which(df3$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),], Tissue) %>% summarize(m = min(OR))
group_by(df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),], Tissue) %>% summarize(m = min(OR))

group_by(df3[which(df3$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),],Tissue) %>% summarize(n = n())
group_by(df3[which(df3$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),], Tissue) %>% summarize(n = n())
group_by(df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),], Tissue) %>% summarize(n = n())

trans<-df3[which(df3$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),]
o<-order(trans$Pvalue)
trans<-trans[o,]


cis<-df2[which(df2$cis_snp==paste0("cis only (N=",nvar[3,2]," SNPs)")),]
a<-(group_by(cis, Type) %>% summarize(n = n()))
dim(a) #172
cis<-cis[which(cis$Pvalue<padj),]
b<-(group_by(cis, Type) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
length(unique(b$Type)) #165

a<-data.frame(full_join(a,b,by="Type"))
a$prop<-a[,3]/a[,2]
#tissue specificity
a[which(a$prop<0.5),]
#      Type n.x n.y         minp  min.OR.  max.OR.      prop
#171 ZNF274   7   1 3.293827e-09 2.029056 2.029056 0.1428571


a[which(a$n.y>20),]
#    Type n.x n.y          minp  min.OR.  max.OR. prop
#32  CTCF  99  99  3.868186e-93 1.951201 2.395581    1
#111 Pol2  52  52 6.942521e-110 2.324978 3.123068    1
a[which(a$n.y>10),]
#        Type n.x n.y          minp  min.OR.  max.OR.      prop
#28     c-Myc  21  20  1.968716e-91 2.132622 3.538772 0.9523810
#32      CTCF  99  99  3.868186e-93 1.951201 2.395581 1.0000000
#50      EZH2  14  11  1.879068e-18 1.788996 3.946814 0.7857143
#102     NRSF  12  12  6.555798e-47 1.845819 2.826837 1.0000000
#103     p300  12  11  2.064187e-58 2.018256 2.542165 0.9166667
#111     Pol2  52  52 6.942521e-110 2.324978 3.123068 1.0000000
#112 Pol2-4H8  12  12 4.134356e-107 2.402047 2.851270 1.0000000
#120    Rad21  12  12  6.769099e-65 1.954174 2.445984 1.0000000
#163      YY1  13  13  7.626693e-74 2.277942 3.074336 1.0000000

b<-(group_by(cis, Type,Tissue) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
dim(b) #391

#blood specific signals
b2<-(group_by(b, Type) %>% summarize(n = n())) #156 TFBS in one tissue
w<-which(b2$n==1)
b2<-data.frame(b[b$Type%in%b2$Type[w],])
b2<-b2[which(b2$Tissue%in%c("Blood")),]


####
amb<-df2[which(df2$cis_snp==paste0("ambivalent (N=",nvar[2,2]," SNPs)")),]
a<-(group_by(amb, Type) %>% summarize(n = n()))
dim(a) #172
amb<-amb[which(amb$Pvalue<padj),]
b<-(group_by(amb, Type) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
length(unique(b$Type)) #156

a<-data.frame(full_join(a,b,by="Type"))
a$prop<-a[,3]/a[,2]
#tissue specificity
a[which(a$prop<0.5),]
#  Type n.x n.y      prop
#50 EZH2  14   6 0.4285714
#62   GR   5   2 0.4000000

a[which(a$n.y>20),]
#     Type n.x n.y          minp  min.OR.  max.OR.     prop
#28  c-Myc  21  21  2.795373e-71 2.152482 3.113537 1.000000
#32   CTCF  99  98  2.127565e-66 1.831038 2.271850 0.989899
#111  Pol2  52  52 1.078096e-110 1.988515 3.065041 1.000000

b<-(group_by(amb, Type,Tissue) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
dim(b) #369

#blood specific signals
b2<-(group_by(b, Type) %>% summarize(n = n())) #156 TFBS in one tissue
w<-which(b2$n==1)
b2<-data.frame(b[b$Type%in%b2$Type[w],])
b2<-b2[which(b2$Tissue%in%c("Blood")),]

#
trans<-df2[which(df2$cis_snp==paste0("trans only (N=",nvar[4,2]," SNPs)")),]
a<-(group_by(trans, Type) %>% summarize(n = n()))
dim(a) #172
trans<-trans[which(trans$Pvalue<padj),]
b<-(group_by(trans, Type) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
length(unique(b$Type)) #156

a<-data.frame(full_join(a,b,by="Type"))
a$prop<-a[,3]/a[,2]
#tissue specificity
a[which(a$prop<0.5),]
#      Type n.x n.y         minp  min.OR.  max.OR.      prop
#171 ZNF274   7   2 4.734275e-05 12.53845 12.94199 0.2857143

a[which(a$n.y>20),]

b<-(group_by(trans, Type,Tissue) %>% summarize(n = n(),minp=min(Pvalue),min(OR),max(OR)))
dim(b) #2

#blood specific signals
b2<-(group_by(b, Type) %>% summarize(n = n())) #156 TFBS in one tissue
w<-which(b2$n==1)
b2<-data.frame(b[b$Type%in%b2$Type[w],])
b2<-b2[which(b2$Tissue%in%c("Blood")),]
