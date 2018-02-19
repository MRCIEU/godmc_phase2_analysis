#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)

filter<-as.numeric(10)

padj<-read.table(paste0(path,"mqtl_segmentations1/garfield.Meff.mqtl_segmentations1.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

pval_lim<-padj

tiss<-read.table("~/repo/godmc_phase2_analysis/07_enrichments/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
ann<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T,sep="\t")

snptype<-c("mqtl_segmentations","mqtl_ambivalent_segmentations","mqtl_cis_segmentations","mqtl_trans_segmentations")
snptype2<-c("All","ambivalent","cis_only","trans_only")

r.all<-data.frame()
for (j in 1:length(snptype)) {
cat(snptype[j],"\n")
for (i in 1:25){
p<-paste0(path,snptype[j],i,"/","garfield.test.",snptype[j],i,".out")
f<-file.size(p)
if(!is.na(f)){
cat(i,"\n")
r<-read.table(paste0(path,snptype[j],i,"/","garfield.test.",snptype[j],i,".out"),he=T)
r$STATE<-ann$NO.[i]
r$snp_cis<-snptype2[j]
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

nvar<-r.all2%>%group_by(snp_cis)%>%summarize(mean(NThresh))
r.all2$snp_cis<-gsub("All",paste0("All (N=",nvar[1,2]," SNPs)"),r.all2$snp_cis)
r.all2$snp_cis<-gsub("ambivalent",paste0("ambivalent (N=",nvar[2,2]," SNPs)"),r.all2$snp_cis)
r.all2$snp_cis<-gsub("cis_only",paste0("cis only (N=",nvar[3,2]," SNPs)"),r.all2$snp_cis)
r.all2$snp_cis<-gsub("trans_only",paste0("trans only (N=",nvar[4,2]," SNPs)"),r.all2$snp_cis)



m<-max(-log10(r.all2$Pvalue))
r.all2$STATE<-gsub("prime","'",r.all2$STATE)

max(r.all2[which(r.all2$STATE%in%c("TxWk","Tx3'","Tx5'")&r.all2$snp_cis%in%c("cis_only","ambivalent")&r.all2$Pvalue<padj),"OR"])

print("repressed")
print(table(r.all2[which(r.all2$OR<1&r.all2$Pvalue<padj),"STATE"]))
print("activated")
print(table(r.all2[which(r.all2$OR>1&r.all2$Pvalue<padj),"STATE"]))


r.all3<-r.all2[which(r.all2$snp_cis!=paste0("All (N=",nvar[1,2]," SNPs)")),]
  
  p1<-ggplot(r.all3,aes(x=STATE,y=OR))+
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
  ggsave(p1,file=paste0("./images/snp_segmentations_OR.pdf"),height=10,width=18)

p1<-ggplot(r.all3,aes(x=STATE,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~snp_cis,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),legend.position="bottom")+
  scale_size(range=c(1,4))+
  ylim(0,(m+0.2*m))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  guides(fill = guide_legend(ncol=10))
ggsave(p1,file=paste0("./images/snp_segmentations.pdf"),height=10,width=18)



#STATE NO. MNEMONIC  DESCRIPTION COLOR NAME  COLOR CODE
#1 TssA  Active TSS  Red 255,0,0
#2 PromU Promoter Upstream TSS Orange Red  255,69,0 ambivalent,cis
#3 PromD1  Promoter Downstream TSS with DNase  Orange Red  255,69,0 ambivalent
#4 PromD2  Promoter Downstream TSS Orange Red  255,69,0
#5 Tx5'  Transcription 5'  Green 0,128,0 ambivalent,trans
#6 Tx  Transcription Green 0,128,0 ambivalent, 
#7 Tx3'  Transcription 3'  Green 0,128,0 trans
#8 TxWk  Weak transcription  Light Green 0,150,0 cis, ambivalent, trans
#39 TxReg Transcription Regulatory  GreenYellow 194,225,5
#10  TxEnh5' Transcription 5' Enhancer GreenYellow 194,225,5
#11  TxEnh3' Transcription 3' Enhancer GreenYellow 194,225,5
#12  TxEnhW  Transcription Weak Enhancer GreenYellow 194,225,5 ambivalent
#13  EnhA1 Active Enhancer 1 Orange  255,195,77
#14  EnhA2 Active Enhancer 2 Orange  255,195,77
#15  EnhAF Active Enhancer Flank Orange  255,195,77
#16  EnhW1 Weak Enhancer 1 Yellow  255,255,0
#17  EnhW2 Weak Enhancer 2 Yellow  255,255,0: ambivalent,cis
#18  EnhAc Enhancer Acetylation Only Yellow  255,255,0
#19  DNase DNase only  Light Yellow  255,255,102
#20  ZNF/Rpts  ZNF genes & repeats Medium Aquamarine 102,205,170
#21  Het Heterochromatin PaleTurquoise 138,145,208
#22  PromP Poised Promoter Light Purple  230,184,183
#23  PromBiv Bivalent Promoter Purple  112,48,160
#24  ReprPC  Repressed PolyComb  Silver  128,128,128 cis
#25  Quies Quiescent/Low White 255,255,255

