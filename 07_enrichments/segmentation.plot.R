#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)

filter<-as.numeric(10)
padj<-read.table(paste0(path,"mqtl_segmentations1/garfield.Meff.mqtl_segmentations1.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

pval_lim<-padj

ann<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T,sep="\t")

snptype<-c("mqtl_segmentations","mqtl_ambivalent_segmentations","mqtl_cis_segmentations","mqtl_trans_segmentations")

for (j in 1:length(snptype)) {
cat(snptype[j],"\n")
r.all<-data.frame()
for (i in 1:25){
p<-paste0(path,snptype[j],i,"/","garfield.test.",snptype[j],i,".out")
f<-file.size(p)
if(!is.na(f)){
cat(i,"\n")
r<-read.table(paste0(path,snptype[j],i,"/","garfield.test.",snptype[j],i,".out"),he=T)
r$STATE<-ann$NO.[i]
r.all<-rbind(r.all,r)
}}
r.all <- unique(r.all[which(as.numeric(as.character(r.all$NThresh))>=filter),])
  r.all$Annotation<-gsub("/lustre/scratch113/projects/uk10k/users/vi1/segmentations/data/NIH_25Marks/","",r.all$Annotation)
  r.all$Annotation<-gsub("_25_imputed12marks_stateno.bed.gz","",r.all$Annotation)

w<-which(is.na(r.all$STATE))
length(w)

tiss<-read.table("~/repo/godmc_phase2_analysis/07_enrichments/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
    m<-match(r.all$Annotation,tiss$ID)
    r.all$Tissue<-tiss[m,"Tissue"]
    r.all$Annotation<-tiss[m,"Name"]
    r.all$Type<-tiss[m,"Type"]
    r.all$Celltype<-tiss[m,"Tissue"]

#

r.all2<-r.all[which(r.all$PThres==1e-14),]
r.all2$logOddsRatio<-log(r.all2$OR)

w<-which(r.all2$Pvalue==0)
m<-min(r.all2$Pvalue[-w])
r.all2$Pvalue[w]<-m
m<-max(-log10(r.all2$Pvalue))
r.all2$STATE<-gsub("prime","'",r.all2$STATE)
p1<-ggplot(r.all2,aes(x=STATE,y=-log10(Pvalue),size=logOddsRatio,fill=Tissue))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),legend.position="bottom")+
  scale_size(range=c(1,4))+
  ylim(0,(m+0.2*m))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))

ggsave(p1,file=paste0(snptype[j],".pdf"),height=10,width=18)

print("repressed")
print(table(r.all2[which(r.all2$OR<1&r.all2$Pvalue<padj),"STATE"]))
print("activated")
print(table(r.all2[which(r.all2$OR>1&r.all2$Pvalue<padj),"STATE"]))

}


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





