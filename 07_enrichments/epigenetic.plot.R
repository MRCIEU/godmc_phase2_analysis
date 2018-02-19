#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)

r<-read.table(paste0(path,"mqtl_epigenetic/garfield.test.mqtl_epigenetic.out"),he=T)
r$cis_snp<-"All"

cis<-read.table(paste0(path,"mqtl_cis_epigenetic/garfield.test.mqtl_cis_epigenetic.out"),he=T)
cis$cis_snp<-"cis only"

trans<-read.table(paste0(path,"mqtl_trans_epigenetic/garfield.test.mqtl_trans_epigenetic.out"),he=T)
trans$cis_snp<-"trans only"

amb<-read.table(paste0(path,"mqtl_ambivalent_epigenetic/garfield.test.mqtl_ambivalent_epigenetic.out"),he=T)
amb$cis_snp<-"ambivalent"

padj<-read.table(paste0(path,"mqtl_epigenetic/garfield.Meff.mqtl_epigenetic.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

df<-rbind(r,cis,trans,amb)
df$logOddsRatio<-log(df$OR)
df2<-df[df$PThresh=="1e-14",]
df2$Type<-as.factor(df2$Type)
df2$Category<-as.factor(df2$Category)
cats<-levels(df2$Category)

pval_lim<-padj

w<-which(df2$Pvalue==0)
m<-min(df2$Pvalue[-w])
df2$Pvalue[w]<-m

for (i in 1:(length(cats))){
cat(i,"\n")

if(cats[i]%in%c("FAIRE","Footprints","Hotspots","Peaks","TFBS")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Tissue,y=-log10(Pvalue),size=logOddsRatio,fill=cis_snp))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
#facet_wrap(~cis_snp,scale="free_x",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/epigenetic_",cats[i],".pdf"),height=10,width=12)
}

if(cats[i]%in%c("Chromatin_States","Histone_Modifications")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
facet_wrap(~cis_snp,scale="free_x",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/epigenetic_",cats[i],".pdf"),height=10,width=12)
}


if(cats[i]%in%c("Genic")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=cis_snp))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  ylim(0,(m+3))+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=15),legend.position="bottom")+
  scale_size(range=c(1,8))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/epigenetic_",cats[i],".pdf"),height=10,width=18)
}
cat(cats[i],"\n")
cat(max(df3[which(df3$cis_snp%in%"All"),"Pvalue"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"trans only"),"Pvalue"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"cis only"),"Pvalue"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"ambivalent"),"Pvalue"]),"\n")

cat(min(df3[which(df3$cis_snp%in%"All"),"OR"]),"\n")
cat(min(df3[which(df3$cis_snp%in%"trans only"),"OR"]),"\n")
cat(min(df3[which(df3$cis_snp%in%"cis only"),"OR"]),"\n")
cat(min(df3[which(df3$cis_snp%in%"ambivalent"),"OR"]),"\n")

cat(max(df3[which(df3$cis_snp%in%"All"),"OR"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"trans only"),"OR"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"cis only"),"OR"]),"\n")
cat(max(df3[which(df3$cis_snp%in%"ambivalent"),"OR"]),"\n")

}

group_by(df3[which(df3$Category=="Peaks"),], Tissue) %>% summarize(m = mean(OR))
group_by(df3[which(df3$Category=="Peaks"),], Tissue) %>% summarize(m = min(OR))
group_by(df3[which(df3$Category=="Peaks"),], Tissue) %>% summarize(m = max(OR))
group_by(df3[which(df3$Category=="Peaks"&df3$Tissue=="fetal_brain"),], Tissue) %>% summarize(m = min(OR))
group_by(df3[which(df3$Category=="Peaks"&df3$Tissue=="fetal_brain"),], Tissue) %>% summarize(m = max(OR))

i=4
df3<-df2[which(df2$Category==cats[i]),c("Annotation","cis_snp","Type","OR","Pvalue","Beta")]df3[which(df3$Pvalue<padj),]
df3[which(df3$Pvalue<padj),]

i=7
df3<-df2[which(df2$Category==cats[i]),c("Annotation","cis_snp","Type","OR","Pvalue","Beta","Tissue")]
table(df3$cis_snp)

#      All ambivalent   cis only trans only 
#      424        424        424        424 
df3<-df3[which(df3$Pvalue<padj),]

table(df3$cis_snp)

#       All ambivalent   cis only trans only 
#       424        424        424        134 
group_by(df3[which(df3$cis_snp=="ambivalent"),]) %>% summarize(m = max(Pvalue))
group_by(df3[which(df3$cis_snp=="cis only"),]) %>% summarize(m = max(Pvalue))
group_by(df3[which(df3$cis_snp=="trans only"),]) %>% summarize(m = max(Pvalue))

group_by(df3[which(df3$cis_snp=="All"),],Tissue) %>% summarize(m = max(OR)) # blood 2.368657
group_by(df3[which(df3$cis_snp=="cis only"),], Tissue) %>% summarize(m = max(OR))#blood 2.227886
group_by(df3[which(df3$cis_snp=="ambivalent"),], Tissue) %>% summarize(m = max(OR))#blood 2.000729

data.frame(group_by(df3[which(df3$cis_snp=="trans only"),], Tissue) %>% summarize(m = min(OR))) #blood 0.4539215

group_by(df3[which(df3$Tissue=="blood"&df3$cis_snp=="ambivalent" ),]) %>% summarize(m = min(OR))#1.654409
group_by(df3[which(df3$Tissue=="blood"&df3$cis_snp=="cis only" ),]) %>% summarize(m = min(OR))#1.889232
group_by(df3[which(df3$Tissue=="blood"&df3$cis_snp=="ambivalent" ),]) %>% summarize(m = max(OR))#2.000729
group_by(df3[which(df3$Tissue=="blood"&df3$cis_snp=="cis only" ),]) %>% summarize(m = max(OR))#2.227886


group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="ambivalent" ),]) %>% summarize(m = min(OR))#1.407696
group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="cis only" ),]) %>% summarize(m = min(OR))#1.594784
group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="trans only" ),]) %>% summarize(m = min(OR))#1.594784

group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="ambivalent" ),]) %>% summarize(m = max(OR))#1.490781
group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="cis only" ),]) %>% summarize(m = max(OR))#1.693076
group_by(df3[which(df3$Tissue=="fetal_brain"&df3$cis_snp=="trans only" ),]) %>% summarize(m = max(OR))#1.594784






