#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)


cis<-read.table(paste0(path,"mqtl_cis_tfbs_sds/garfield.test.mqtl_cis_tfbs_sds.out"),he=T)
#cis$snp_cis<-"ambivalent+cis only"

cis2<-read.table(paste0(path,"cismqtl_sds_extremes/garfield.test.cismqtl_sds_extremes.out"),he=T)
#cis2$snp_cis<-"ambivalent+cis only"
cis2$Category<-gsub("TFBS","TFBS_celltype",cis2$Category)

padj<-read.table(paste0(path,"mqtl_cis_tfbs_sds/garfield.Meff.mqtl_cis_tfbs_sds.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

df<-rbind(cis,cis2)

#df2%>%group_by(snp_cis)%>%summarize(mean(NThresh))
df$logOddsRatio<-log(df$OR)
df2<-df[df$PThresh=="1e-14",]
df2$snp_cis<-paste0("ambivalent+cis only (N=",round(mean(df2$NThresh),0)," SNPs)")
df2$Type<-as.factor(df2$Type)
df2$Category<-as.factor(df2$Category)
cats<-levels(df2$Category)

pval_lim<-padj

w<-which(df2$Pvalue==0)
if(length(w)>0){
m<-min(df2$Pvalue[-w])
df2$Pvalue[w]<-m}

for (i in 1:(length(cats))){
cat(i,"\n")

if(cats[i]%in%c("TFBS")) {
df3<-df2[which(df2$Category==cats[i]),]
df3$Type<-sub("_[^_]+$", "", df3$Type)
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
facet_wrap(~snp_cis,scale="free_x",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/sds_",cats[i],".pdf"),height=10,width=16)
}


if(cats[i]%in%c("FAIRE","Footprints","Hotspots","Peaks","TFBS_celltype")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Tissue,y=-log10(Pvalue),size=logOddsRatio,fill=snp_cis))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
facet_wrap(~snp_cis,scale="free_x",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/sds_",cats[i],".pdf"),height=10,width=12)
}

if(cats[i]%in%c("Chromatin_States","Histone_Modifications")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
geom_point(alpha=0.7,shape=21,stroke=1)+
ylim(0,(m+50))+
facet_wrap(~snp_cis,scale="free_x",ncol=1)+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
scale_size(range=c(1,4))+
scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/sds_",cats[i],".pdf"),height=10,width=12)
}


if(cats[i]%in%c("Genic")) {
df3<-df2[which(df2$Category==cats[i]),]
m<-max(-log10(df3$Pvalue))
p1<-ggplot(df3,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=snp_cis))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  ylim(0,(m+3))+
  facet_wrap(~snp_cis,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=15),legend.position="bottom")+
  scale_size(range=c(1,8))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
ggsave(p1,file=paste0("./images/sds_",cats[i],".pdf"),height=10,width=18)
}
}

##
tfbs<-df2[which(df2$Category=="TFBS"),]
tfbs$Type<-sub("_[^_]+$", "", tfbs$Type)
tfbs<-tfbs[which(tfbs$Pvalue<padj),]
o<-order(tfbs$Pvalue)
tfbs<-tfbs[o,]



filter<-as.numeric(10)

padj<-read.table(paste0(path,"mqtl_cis_sds_segmentations1/garfield.Meff.mqtl_cis_sds_segmentations1.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

pval_lim<-padj

tiss<-read.table("~/repo/godmc_phase2_analysis/07_enrichments/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
ann<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T,sep="\t")

#snptype<-c("mqtl_segmentations","mqtl_ambivalent_segmentations","mqtl_cis_segmentations","mqtl_trans_segmentations")
#snptype2<-c("All","ambivalent","cis_only","trans_only")


r.all<-data.frame()
#for (j in 1:length(snptype)) {
#cat(snptype[j],"\n")
for (i in 1:25){
p<-paste0(path,"mqtl_cis_sds_segmentations",i,"/","garfield.test.mqtl_cis_sds_segmentations",i,".out")
f<-file.size(p)
if(!is.na(f)){
cat(i,"\n")
r<-read.table(paste0(path,"mqtl_cis_sds_segmentations",i,"/","garfield.test.mqtl_cis_sds_segmentations",i,".out"),he=T)
r$STATE<-ann$NO.[i]
r$snp_cis<-paste0("ambivalent+cis only (N=",round(mean(r$NThresh),0)," SNPs)")
r$Annotation<-gsub("/lustre/scratch113/projects/uk10k/users/vi1/segmentations/data/NIH_25Marks/","",r$Annotation)
r$Annotation<-gsub("_25_imputed12marks_stateno.bed.gz","",r$Annotation)

m<-match(r$Annotation,tiss$ID)
    r$Tissue<-tiss[m,"Tissue"]
    r$Annotation<-tiss[m,"Name"]
    r$Type<-tiss[m,"Type"]
    r$Celltype<-tiss[m,"Tissue"]
r.all<-rbind(r.all,r)


}}
#r.all <- unique(r.all[which(as.numeric(as.character(r.all$NThresh))>=filter),])
  
w<-which(is.na(r.all$STATE))
length(w)

    
#

r.all2<-r.all[which(r.all$PThres==1e-14),]
r.all2$logOddsRatio<-log(r.all2$OR)

w<-which(r.all2$Pvalue==0)
if (length(w)>0){
m<-min(r.all2$Pvalue[-w])
r.all2$Pvalue[w]<-m
}
}

m<-max(-log10(r.all2$Pvalue))
r.all2$STATE<-gsub("prime","'",r.all2$STATE)


print("repressed")
print(table(r.all2[which(r.all2$OR<1&r.all2$Pvalue<padj),"STATE"]))
print("activated")
print(table(r.all2[which(r.all2$OR>1&r.all2$Pvalue<padj),"STATE"]))


#max(r.all2[which(r.all2$STATE%in%c("TxWk","Tx3'","Tx5'")&r.all2$snp_cis%in%c("cis_only","ambivalent")&r.all2$Pvalue<padj),"OR"])

p1<-ggplot(r.all2,aes(x=STATE,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~snp_cis,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),legend.position="bottom")+
  scale_size(range=c(1,4))+
  ylim(0,(m+0.2*m))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  guides(fill = guide_legend(ncol=10))
ggsave(p1,file=paste0("./images/snp_sds_segmentations.pdf"),height=10,width=18)


  
  p1<-ggplot(r.all2,aes(x=STATE,y=OR))+
  geom_hline(yintercept=1, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=Tissue,size=-log10(Pvalue)))+
  facet_wrap(~snp_cis,ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  #pvalue scaling
  scale_size(range=c(1,8))+
  guides(fill = guide_legend(ncol=20))+
  scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,6),limits=c(0.2,6)) +
  ylab("Odds ratio (log scale)") +
  theme(legend.text=element_text(size=12))
  ggsave(p1,file=paste0("./images/snp_sds_segmentations_OR.pdf"),height=10,width=18)


###


