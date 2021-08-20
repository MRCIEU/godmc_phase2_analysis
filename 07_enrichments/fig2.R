#module add languages/R-3.5.1-ATLAS-gcc-6.1
path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(scales)

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
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
r.all2<-r.all[which(r.all$PThres==1e-14),]
r.all2$logOR<-log(r.all2$OR)
r.all2$abslogOR <- abs(log(r.all2$OR))
r.all2$snp_cis<-gsub("ambivalent","cis+trans",r.all2$snp_cis)
r.all2$snp_cis<-gsub("cis_only","cis only",r.all2$snp_cis)
r.all2$snp_cis<-gsub("trans_only","trans only",r.all2$snp_cis)
r.all2$FDRP<-p.adjust(r.all2$Pvalue,method="BH")
w<-which(r.all2$snp_cis%in%c("All"))

p1<-ggplot(r.all2[-w,], aes(x=logOR)) +
  geom_density(aes(fill=snp_cis), alpha=0.5, bw=0.1) +
  labs(fill = "SNP annotation") +
  scale_fill_brewer(type="qual") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 3)

ggsave(p1,file="fig2_segmentation_snps.pdf")

bl<-r.all2[which(r.all2$Tissue=="BLOOD"&r.all2$Type=="PrimaryCell"),]
g<-grep("cord",bl$Annotation)
bl<-bl[-g,]
w<-which(bl$snp_cis%in%c("All"))

ann<-data.frame(table(bl[w,"Annotation"]))
table(ann$Freq)
#0 25
#105 22
###
bl$STATE<-as.character(bl$STATE)

snp<-bl%>%group_by(STATE,snp_cis)%>%summarize(OR=mean(OR),Pval=mean(Pvalue),FDRp=mean(FDRP))
snp1<-data.frame(feature="SNP",snp)
names(snp1)<-c("feature","STATE","Annotation","OR","Pvalue","FDR Pvalue")
snp<-spread(snp,key=snp_cis,value=OR)

#snp<-reshape(snp[,c("STATE","snp_cis","OR")],idvar="STATE",timevar = "snp_cis", direction = "wide")

###

load("../results/enrichments/lola_chromstates_mqtlcpg.rdata")
tiss<-read.table("jul2013.roadmapData_tissues.txt",sep="\t",he=T)


LOLA_res<-lola_res0_matched
LOLA_res$userSet<-as.character(LOLA_res$userSet)
LOLA_res$userSet<-"All"

#load("../results/enrichments/lola_chromstates_mqtlcpg_cis.rdata")
load("../results/enrichments/lola_chromstates_mqtlcpg_cis_updated.rdata")
w<-which(res_all$userSet=="ambivalent (10941 regions)")
res_all$userSet[w]<-"cis+trans"
w<-which(res_all$userSet=="trans_only (3864 regions)")
res_all$userSet[w]<-"trans only"
w<-which(res_all$userSet=="cis_only (107004 regions)")
res_all$userSet[w]<-"cis only"

LOLA_res<-rbind(LOLA_res,res_all)
LOLA_res<-LOLA_res[which(LOLA_res$collection%in%c("25states")),]
LOLA_res$logOR <- log(LOLA_res$oddsRatio)
LOLA_res$abslogOR <- abs(log(LOLA_res$oddsRatio))
LOLA_res<-LOLA_res[which(LOLA_res$logOddsRatio!="-Inf"),]
LOLA_res[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
LOLA_res[,Pvalue:=10^(-pValueLog)]

m<-match(LOLA_res$dataSource,tiss$ID)
LOLA_res<-data.frame(LOLA_res,tiss[m,])


w<-which(LOLA_res$userSet%in%c("All"))

p2<-ggplot(LOLA_res[-w,], aes(x=logOR)) +
  geom_density(aes(fill=userSet), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  labs(fill = "CpG annotation") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 3)

ggsave(p2,file="fig2_segmentation_cpgs_updated.pdf",device="pdf",height=6,width=6)

LOLA_res$blood<-as.character(LOLA_res$tissue)
w<-which(LOLA_res$blood%in%c("BLOOD")==F)
LOLA_res$blood[w]<-"Other"
LOLA_res$blood<-gsub("BLOOD","Blood",LOLA_res$blood)

LOLA_res$blood2<-LOLA_res$blood
w<-which(LOLA_res$tissue%in%c("BRAIN"))
LOLA_res$blood2[w]<-"Brain"
w<-which(LOLA_res$tissue%in%c("FAT"))
LOLA_res$blood2[w]<-"Fat"

w<-which(LOLA_res$userSet%in%c("All"))
p2<-ggplot(LOLA_res[-w,], aes(x=logOR)) +
  geom_density(aes(fill=blood), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 3)

ggsave(p2,file="fig2_segmentation_cpgs_updated_tissue.pdf")

w<-which(LOLA_res$userSet%in%c("All"))
p2<-ggplot(LOLA_res[-w,], aes(x=logOR)) +
  geom_density(aes(fill=Group), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 3)

ggsave(p2,file="fig2_segmentation_cpgs_updated_tissuegroup.pdf")


LOLA_res2<-LOLA_res[-w,]

g2<-grep("Enh",LOLA_res2$seg_explanation)
LOLA_enh<-LOLA_res2[g2,]
tissues<-unique(LOLA_enh$tissue)

for (i in 1:length(tissues)){
g1<-which(LOLA_enh$tissue%in%c("BLOOD",as.character(tissues[i])))
p2<-ggplot(LOLA_enh[g1,], aes(x=logOR)) +
  geom_density(aes(fill=tissue),alpha=0.5, bw=0.1) +
  #scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)

ggsave(p2,file=paste0("fig2_segmentation_cpgs_updated_tissues",tissues[i],"_Enh.pdf"),height=6, width=8)
}

g1<-which(LOLA_enh$dataSource%in%c("E003","E011","E012","E013"))
LOLA_es<-LOLA_enh[g1,]

p2<-ggplot(LOLA_es, aes(x=logOR)) +
  geom_density(aes(fill=dataSource),alpha=0.5, bw=0.1) +
  #scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "dataSource") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)

ggsave(p2,file=paste0("fig2_segmentation_cpgs_updated_tissues_es_Enh.pdf"),height=6, width=8)



states<-unique(LOLA_res2$seg_explanation)

LOLA_res3<-LOLA_res2[which(LOLA_res2$blood2%in%c("Blood","Fat","Brain")),]
for (i in 1:length(states)){
g1<-which(LOLA_res2$seg_explanation%in%states[i])
p2<-ggplot(LOLA_res2[g1,], aes(x=logOR)) +
  geom_density(aes(fill=blood), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Blood") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)
lab<-states[i]
if(i==21) {lab<-"ZNF_Rpts"}
ggsave(p2,file=paste0("fig2_segmentation_cpgs_updated_tissue_",lab,".pdf"),height=6, width=8)
}


g1<-grep("Prom",LOLA_res3$seg_explanation)
p2<-ggplot(LOLA_res3[g1,], aes(x=logOR)) +
  geom_density(aes(fill=blood2), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)

ggsave(p2,file="fig2_segmentation_cpgs_updated_fatbrain_Prom.pdf",height=6, width=8)

g2<-grep("Enh",LOLA_res3$seg_explanation)
p2<-ggplot(LOLA_res3[g2,], aes(x=logOR)) +
  geom_density(aes(fill=blood2), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)

ggsave(p2,file="fig2_segmentation_cpgs_updated_fatbrain_Enh.pdf",height=6, width=8)

g3<-grep("Het",LOLA_res3$seg_explanation)
p2<-ggplot(LOLA_res3[g3,], aes(x=logOR)) +
  geom_density(aes(fill=blood2), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  facet_wrap(~userSet,nrow=1) +
  labs(fill = "Tissue") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-2, 2)

ggsave(p2,file="fig2_segmentation_cpgs_updated_fatbrain_Het.pdf",height=6, width=8)



w<-which(LOLA_res$userSet%in%c("cis only"))

chromstate<-LOLA_res[w,] %>%
    group_by(seg_explanation,blood) %>%
    dplyr::summarize(meanOR = mean(oddsRatio, na.rm=TRUE))

chromstate_b<-chromstate[which(chromstate$blood=="Blood"),]
chromstate_o<-chromstate[which(chromstate$blood=="Other"),]

ratio<-(chromstate_b$meanOR/chromstate_o$meanOR)
w<-which(ratio>1.1)
data.frame(chromstate_b[w,],chromstate_o[w,],ratio[w])


bl<-LOLA_res[which(LOLA_res$Tissue=="BLOOD"&LOLA_res$Type=="PrimaryCell"),]
g<-grep("cord",bl$Name)
bl<-bl[-g,]

w<-which(bl$userSet%in%c("All"))
ann<-data.frame(table(bl[w,"Name"]))
table(ann$Freq)
#0 25
#105 22
###
bl$seg_explanation<-as.character(bl$seg_explanation)
cpg1<-bl[,c("userSet","oddsRatio","Pvalue","p.adjust")]

cpg<-bl%>%group_by(seg_explanation,userSet)%>%summarize(OR=mean(oddsRatio),Pval=mean(Pvalue), FDRPval=mean(p.adjust))
cpg1<-data.frame(feature="DNAm site",cpg)
names(cpg1)<-c("feature","STATE","Annotation","OR","Pvalue","FDR Pvalue")

cpg<-spread(cpg,key=userSet,value=OR)
names(cpg)[1]<-"STATE"
comb<-merge(cpg,snp,by="STATE")


comb1<-rbind(cpg1,snp1)
comb1$label<-paste(comb1$feature,comb1$Annotation,sep="_")
comb1$STATE<-gsub("prime","'",comb1$STATE)


snp_gene<-read.table("../23_enrichments_geneannotations/TableSXX_garfield_gene_annotation.txt",he=T,sep="\t")
snp_gene<-data.frame(feature="SNP",snp_gene[,c("Gene.Annotation","Annotation","OR","Pvalue")])
snp_gene$p.adjust<-p.adjust(snp_gene$Pvalue)
names(snp_gene)<-c("feature","STATE","Annotation","OR","Pvalue","FDR Pvalue")
snp_gene$STATE<-gsub("3UTRs","3'UTR",snp_gene$STATE)
snp_gene$STATE<-gsub("5UTRs","5'UTR",snp_gene$STATE)


cpg_gene<-read.table("../23_enrichments_geneannotations/TableSXX_LOLA_gene_annotation_updated.txt",he=T,sep="\t")
cpg_gene$Gene.Annotation<-gsub("3UTRs","3'UTR",cpg_gene$Gene.Annotation)
cpg_gene$Gene.Annotation<-gsub("5UTRs","5'UTR",cpg_gene$Gene.Annotation)
cpg_gene$Gene.Annotation<-gsub("intergenic","Intergenic",cpg_gene$Gene.Annotation)
cpg_gene$Gene.Annotation<-gsub("introns","Intron",cpg_gene$Gene.Annotation)
cpg_gene$Gene.Annotation<-gsub("exons","Exon",cpg_gene$Gene.Annotation)
cpg_gene<-data.frame(feature="DNAm site",cpg_gene[,c("Gene.Annotation","Annotation","OR","Pvalue")])
cpg_gene$p.adjust<-p.adjust(cpg_gene$Pvalue)
names(cpg_gene)<-c("feature","STATE","Annotation","OR","Pvalue","FDR Pvalue")

w<-which(snp_gene$STATE%in%cpg_gene$STATE)
snp_gene<-snp_gene[w,]
w<-which(cpg_gene$STATE%in%snp_gene$STATE)
cpg_gene<-cpg_gene[w,]

gene<-rbind(snp_gene,cpg_gene)

gene$label<-paste(gene$feature,gene$Annotation,sep="_")
comb1<-rbind(gene,comb1)
names(comb1)<-gsub("FDR Pvalue","FDR_Pvalue",names(comb1))
comb1$stars <- cut(comb1$FDR_Pvalue, breaks=c(Inf, 1e-50,1e-10,0.001,-Inf), label=c("***", "**", "*",""))  # Create column of significance labels

w<-which(comb1$Annotation%in%c("All"))
p1 <- ggplot(comb1[-w,], aes(Annotation,as.factor(STATE),fill = log(OR))) +
geom_tile(colour="gray80") +
#theme_bw(10) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
labs(x="mQTL annotation",y="Feature") +
geom_text(aes(label=stars), color="black", size=3) +
scale_fill_gradient2(low = muted("darkblue"), mid = "white", high = muted("red"), midpoint = 0, space = "Lab", na.value = "grey10", guide = "colourbar", limits=c(-1.5, 1.5)) +
facet_wrap(~feature,ncol=2) +
coord_equal()
ggsave(p1,file="chrom_states_heatmap.pdf")
#scale_fill_gradient2(low = muted("darkblue"), mid = "white", high = muted("red"), midpoint = 0, space = "Lab", na.value = "grey10", guide = "colourbar", limits=c(-0.049, 0.049)) +
#
save(comb1,file="chromstates_fig.Robj")
comb2<-merge(snp,cpg,by.x="STATE",by.y="STATE")

comb2a<-spread(comb1,key=label,value=OR)
comb2b<-spread(comb1,key=label,value=stars)



#scale_fill_gradient2(low = muted("darkblue"), mid = "white", high = muted("red"), midpoint = 0, space = "Lab", na.value = "grey10", guide = "colourbar", limits=c(-0.049, 0.049)) +
#


####
####
####
r<-read.table(paste0(path,"mqtl_tfbs/garfield.test.mqtl_tfbs.out"),he=T)
r$cis_snp<-"All"
length(unique(r$Type))
cis<-read.table(paste0(path,"mqtl_cis_tfbs/garfield.test.mqtl_cis_tfbs.out"),he=T)
cis$cis_snp<-"cis only"
length(unique(cis$Type))
trans<-read.table(paste0(path,"mqtl_trans_tfbs/garfield.test.mqtl_trans_tfbs.out"),he=T)
trans$cis_snp<-"trans only"
amb<-read.table(paste0(path,"mqtl_ambivalent_tfbs/garfield.test.mqtl_ambivalent_tfbs.out"),he=T)
amb$cis_snp<-"cis+trans"

r<-rbind(r,cis,amb,trans)
r<-r[r$PThresh=="1e-14",]
r$Type<-gsub("_.*","",r$Type)
length(unique(r$Type)) #171

r$logOR<-log(r$OR)
r$abslogOR <- abs(log(r$OR))
w<-which(r$cis_snp%in%c("All"))

df1<-r[-w,c("logOR","cis_snp")]

p3<-ggplot(r[-w,], aes(x=logOR)) +
geom_density(aes(fill=cis_snp), alpha=0.5, bw=0.1) +
labs(fill = "SNP annotation") +
scale_fill_brewer(type="qual") +
geom_vline(xintercept=0, linetype="dotted") +
theme(legend.position="bottom") +
labs(x="Enrichment (log OR)") +
theme_bw() +
theme(text = element_text(size=12)) +
xlim(-3, 2.75) +
ylim(0, 3)

ggsave(plot=p3,file="fig2_tfbs_snps.pdf")


load("../results/enrichments/lola_core_mqtlcpg.rdata")
LOLA_res<-lola_res0_matched
LOLA_res$userSet<-as.character(LOLA_res$userSet)
LOLA_res$userSet<-"All"

load("../results/enrichments/lola_core_mqtlcpg_cis_updated.rdata")

w<-which(res_all$userSet=="ambivalent (10941 regions)")
res_all$userSet[w]<-"cis+trans"
w<-which(res_all$userSet=="trans_only (3864 regions)")
res_all$userSet[w]<-"trans only"
w<-which(res_all$userSet=="cis_only (107004 regions)")
res_all$userSet[w]<-"cis only"

LOLA_res<-rbind(LOLA_res,res_all)
LOLA_res<-LOLA_res[which(LOLA_res$collection%in%c("encode_tfbs","codex")),]
  
LOLA_res$logOR <- log(LOLA_res$oddsRatio)
LOLA_res$abslogOR <- abs(log(LOLA_res$oddsRatio))
length(unique(LOLA_res$antibody))
#262
LOLA_res$antibody<-gsub("_.*","",LOLA_res$antibody)
length(unique(LOLA_res$antibody))
#228
w<-which(LOLA_res$userSet%in%c("All"))

df<-LOLA_res[-w,c("userSet","logOR")]

p4<-ggplot(LOLA_res[-w,], aes(x=logOR)) +
geom_density(aes(fill=userSet), alpha=0.5, bw=0.1) +
scale_fill_brewer(type="qual") +
labs(fill = "CpG annotation") +
geom_vline(xintercept=0, linetype="dotted") +
theme(legend.position="bottom") +
labs(x="Enrichment (log OR)") +
theme_bw() +
theme(text = element_text(size=12)) +
xlim(-3, 2.75) +
ylim(0, 3)

#ggsave(p4,file="fig2_tfbs_cpgs.pdf")
ggsave(p4,file="fig2_tfbs_cpgs_updated.pdf")

pdf("Fig2b_tfbs_updated.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p3, vp = vplayout(1, 1))
print(p4, vp = vplayout(1, 2))

dev.off()

df1<-data.frame(ciscat=df1$cis_snp,logOR=df1$logOR,annotation="TFBS amongst mQTL SNPs")
df<-data.frame(ciscat=df$userSet,logOR=df$logOR,annotation="TFBS amongst mQTL CpGs")
df2<-rbind(df,df1)

p5<-ggplot(df2, aes(x=logOR)) +
  geom_density(aes(fill=ciscat), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  labs(fill = "Annotation") +
  facet_wrap(~annotation,nrow=1)+
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  theme_bw() +
  xlim(-3, 2.75) +
  ylim(0, 3)
ggsave(p5,file="fig2_tfbs_facet_updated.pdf",dpi=1200)

###

####
#gwa<-read.csv("GARFIELD_CisAmbvTrans.csv",he=T) #1095=219 traits
gwa<-read.table("gwa_traits.txt",he=T,sep="\t")#462
names(gwa)<-gsub("mQTL.enrichment.category","Annotation",names(gwa))
gwa$Annotation<-gsub("ambv","cis+trans",gwa$Annotation)
gwa$Annotation<-gsub("cis","cis only",gwa$Annotation)
gwa$Annotation<-gsub("trans","trans only",gwa$Annotation)
gwa$Annotation<-gsub("cis only+trans only","cis+trans",gwa$Annotation,fixed=T)

w<-which(gwa$Annotation%in%c("cis only","cis+trans","trans only"))
gwa<-gwa[w,] #462
#w<-which(is.na(gwa$OR))
#gwa$OR[w]<-1e-10
gwa$logOR <- log(gwa$OR)
gwa$abslogOR <- abs(log(gwa$OR))

g<-grep("Kettunen",gwa$Label)#339
gwa<-gwa[-g,]#123
length(unique(gwa$Phenotype)) #37
length(unique(gwa$Label)) #43

#g1<-grep("HDL.",gwa$gwa)
#g2<-grep("LDL.",gwa$gwa)
#g3<-grep("VLDL.",gwa$gwa)
#g4<-grep("IDL.",gwa$gwa)

#g<-c(g1,g2,g3,g4) #112 traits
#gwa<-gwa[-g,]
#length(unique(gwa$gwa))
#[1] 140

#w<-which(gwa$gwas%in%c("Ace","AcAce","Ala","Alb","bOHBut","Bis.DB.ratio","Bis.FA.ratio","Cit","CH2.DB.ratio","CH2.in.FA","DB.in.FA","DHA","Est.C","Free.C","FAw3","FAw6","FAw67","FAw79S","Glc","Gln","Gly","Gp","His","Ile","LA","Lac","Leu","MUFA","Phe","Pyr","PC","otPUFA","Tyr","Tot.FA","Val","Serum.C","Serum.TG","SM","TotPG"))
#gwa<-gwa[-w,]
#length(unique(gwa$gwa)) #102

#w<-which(gwa$gwas%in%c("SDS"))
#gwa<-gwa[-w,]
#length(unique(gwa$gwa)) #101

#duplicated GWAS
#w<-which(gwa$gwas%in%c("CAD_2014","Alzheimers","ALS_lmm"))
#gwa<-gwa[-w,]

#length(unique(gwa$gwa)) #102

gwa %>%
    group_by(Annotation) %>%
    dplyr::summarize(Mean = mean(logOR, na.rm=TRUE))

# A tibble: 3 x 2
#  Annotation       Mean
#       <chr>      <dbl>
#1   cis only  0.9511548
#2  cis+trans -2.2809449
#3 trans only -8.8027932

# A tibble: 3 x 2
#Annotation        Mean
#<chr>       <dbl>
#  1   cis only   0.6871599
#2  cis+trans   0.7491646
#3 trans only -11.3796267

nthresh<-gwa %>%
    group_by(Phenotype) %>%
    dplyr::summarize(sum = sum(NThresh, na.rm=TRUE))

df<-data.frame(nthresh)
w<-which(df$sum==0)
length(w) #0

#w<-which(gwa$Phenotype%in%df[w,1]) #81
#gwa<-gwa[-w,] #traits
#length(unique(gwa$gwa)) #74

gwa %>%
    group_by(Annotation) %>%
    dplyr::summarize(Mean = mean(logOR, na.rm=TRUE))
# A tibble: 3 x 2
#Annotation        Mean
#<chr>       <dbl>
#  1   cis only   0.6871599
#2  cis+trans   0.7491646
#3 trans only -11.3796267
nthresh<-gwa %>%
    group_by(Phenotype) %>%
    dplyr::summarize(sum = sum(NThresh, na.rm=TRUE))

df<-data.frame(nthresh)
w<-which(df$sum<10)
length(w) #0
#w<-which(gwa$gwas%in%df[w,1])
#gwa<-gwa[-w,]#33
#length(unique(gwa$gwa)) #63

gwa$FDR<-p.adjust(gwa$Pvalue,method="BH")
write.table(gwa,"GARFIELD_gwa_traits.txt",sep="\t",row.names=F,col.names=T,quote=F)

gwa2<-gwa[which(gwa$FDR<0.05),]
table(gwa2$Annotation)

#cis only  cis+trans trans only
#14         19          1
gwa<-data.frame(gwa,gwa="mQTL SNPs amongst GWA SNPs")
p5<-ggplot(gwa, aes(x=logOR)) +
geom_density(aes(fill=Annotation), alpha=0.5, bw=0.1,adjust=3) +
scale_fill_brewer(type="qual") +
geom_vline(xintercept=0, linetype="dotted") +
facet_wrap(~gwa,nrow=1)+
theme(legend.position="bottom") +
theme_bw() +
xlim(-6, 6) +
ylim(0, 1) +
labs(x="Enrichment (log OR)")
ggsave(p5,file="fig2_gwa.pdf",height=6,width=4,dpi=1200)

p5<-ggplot(df2, aes(x=logOR)) +
  geom_density(aes(fill=ciscat), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  labs(fill = "Annotation") +
  facet_wrap(~annotation,nrow=1)+
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 2.75) +
  ylim(0, 3)
ggsave(p5,file="fig2_tfbs_facet.pdf")


p5<-ggplot(gwa, aes(x=abs(logOR))) +
geom_density(aes(fill=Annotation), alpha=0.5, bw=0.1) +
scale_fill_brewer(type="qual") +
geom_vline(xintercept=0, linetype="dotted") +
theme(legend.position="bottom") +
xlim(-3, 3)
ggsave(p5,file="fig2_gwa_absOR.pdf")


sel<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/selection_no_mhc_lct_forplot.txt",he=T)
w<-which(sel$Type%in%c("ambivalentmqtls","cismqtls","transmqtls"))
sel<-sel[w,]
sel$Type<-gsub("ambivalentmqtls","cis+trans",sel$Type)
sel$Type<-gsub("cismqtls","cis only",sel$Type)
sel$Type<-gsub("transmqtls","trans only",sel$Type)

sel<-sel[sel$PThresh=="1e-14",]

sel$logOR<-log(sel$OR)
sel$abslogOR <- abs(log(sel$OR))

p6<-ggplot(sel, aes(x=logOR)) +
geom_density(aes(fill=Type), alpha=0.5, bw=0.1) +
scale_fill_brewer(type="qual") +
geom_vline(xintercept=0, linetype="dotted") +
theme(legend.position="bottom") +
xlim(-3, 3)

ggsave(p6,file="fig2_sel.pdf")

pdf("Fig2.pdf",height=10,width=16)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
print(p3, vp = vplayout(2, 1))
print(p4, vp = vplayout(2, 2))
print(p5, vp = vplayout(3, 1))
dev.off()

###
#module add languages/R-3.5.1-ATLAS-gcc-6.1
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(scales)

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"

r<-read.table(paste0(path,"mqtl_tfbs/garfield.test.mqtl_tfbs.out"),he=T)
r$cis_snp<-"All"
length(unique(r$Type))
cis<-read.table(paste0(path,"mqtl_cis_tfbs/garfield.test.mqtl_cis_tfbs.out"),he=T)
cis$cis_snp<-"cis only"
length(unique(cis$Type))
trans<-read.table(paste0(path,"mqtl_trans_tfbs/garfield.test.mqtl_trans_tfbs.out"),he=T)
trans$cis_snp<-"trans only"
amb<-read.table(paste0(path,"mqtl_ambivalent_tfbs/garfield.test.mqtl_ambivalent_tfbs.out"),he=T)
amb$cis_snp<-"cis+trans"

r<-rbind(r,cis,amb,trans)
r<-r[r$PThresh=="1e-14",]
r$Type<-gsub("_.*","",r$Type)
length(unique(r$Type)) #171

r$logOR<-log(r$OR)
r$abslogOR <- abs(log(r$OR))
w<-which(r$cis_snp%in%c("All"))

df1<-r[-w,c("logOR","cis_snp")]

#p3<-ggplot(r[-w,], aes(x=logOR)) +
#geom_density(aes(fill=cis_snp), alpha=0.5, bw=0.1) +
#labs(fill = "SNP annotation") +
#scale_fill_brewer(type="qual") +
#geom_vline(xintercept=0, linetype="dotted") +
#theme(legend.position="bottom") +
#labs(x="Enrichment (log OR)") +
#xlim(-3, 2.75) +
#ylim(0, 3)

#ggsave(plot=p3,file="fig2_tfbs_snps.pdf")

load("../results/enrichments/lola_core_mqtlcpg.rdata")
LOLA_res<-lola_res0_matched
LOLA_res$userSet<-as.character(LOLA_res$userSet)
LOLA_res$userSet<-"All"

load("../results/enrichments/lola_core_mqtlcpg_cis_updated.rdata")

w<-which(res_all$userSet=="ambivalent (10941 regions)")
res_all$userSet[w]<-"cis+trans"
w<-which(res_all$userSet=="trans_only (3864 regions)")
res_all$userSet[w]<-"trans only"
w<-which(res_all$userSet=="cis_only (107004 regions)")
res_all$userSet[w]<-"cis only"

LOLA_res<-rbind(LOLA_res,res_all)
LOLA_res<-LOLA_res[which(LOLA_res$collection%in%c("encode_tfbs","codex")),]
  
LOLA_res$logOR <- log(LOLA_res$oddsRatio)
LOLA_res$abslogOR <- abs(log(LOLA_res$oddsRatio))
length(unique(LOLA_res$antibody))
#262
LOLA_res$antibody<-gsub("_.*","",LOLA_res$antibody)
length(unique(LOLA_res$antibody))
#228
w<-which(LOLA_res$userSet%in%c("All"))

df<-LOLA_res[-w,c("userSet","logOR")]

#p4<-ggplot(LOLA_res[-w,], aes(x=logOR)) +
#geom_density(aes(fill=userSet), alpha=0.5, bw=0.1) +
#scale_fill_brewer(type="qual") +
#labs(fill = "CpG annotation") +
#geom_vline(xintercept=0, linetype="dotted") +
#theme(legend.position="bottom") +
#labs(x="Enrichment (log OR)") +
#xlim(-3, 2.75) +
#ylim(0, 3)

#ggsave(p4,file="fig2_tfbs_cpgs.pdf")
#ggsave(p4,file="fig2_tfbs_cpgs_updated.pdf")

#pdf("Fig2b_tfbs_updated.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(1, 2)))
#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(p3, vp = vplayout(1, 1))
#print(p4, vp = vplayout(1, 2))

#dev.off()

df1<-data.frame(ciscat=df1$cis_snp,logOR=df1$logOR,annotation="TFBS amongst mQTL SNPs")
df<-data.frame(ciscat=df$userSet,logOR=df$logOR,annotation="TFBS amongst mQTL CpGs")
df2<-rbind(df,df1)

p5<-ggplot(df2, aes(x=logOR)) +
  geom_density(aes(fill=ciscat), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  labs(fill = "Annotation") +
  facet_wrap(~annotation,nrow=1)+
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 2.75) +
  ylim(0, 3)
ggsave(p5,file="fig2_tfbs_facet_updated.pdf")

###

####
#gwa<-read.csv("GARFIELD_CisAmbvTrans.csv",he=T) #1095=219 traits
gwa<-read.table("gwa_traits.txt",he=T,sep="\t")#462
g<-grep("Enigma",gwa$Label)
gwa[g,"Label"]<-"EnigmaPutamen_Hibar_2016"
names(gwa)<-gsub("mQTL.enrichment.category","Annotation",names(gwa))
gwa$Annotation<-gsub("ambv","cis+trans",gwa$Annotation)
gwa$Annotation<-gsub("cis","cis only",gwa$Annotation)
gwa$Annotation<-gsub("trans","trans only",gwa$Annotation)
gwa$Annotation<-gsub("cis only+trans only","cis+trans",gwa$Annotation,fixed=T)

w<-which(gwa$Annotation%in%c("cis only","cis+trans","trans only"))
gwa<-gwa[w,] #462
#w<-which(is.na(gwa$OR))
#gwa$OR[w]<-1e-10
gwa$logOR <- log(gwa$OR)
gwa$abslogOR <- abs(log(gwa$OR))

g<-grep("Kettunen",gwa$Label)#339
gwa<-gwa[-g,]#123
length(unique(gwa$Phenotype)) #37
length(unique(gwa$Label)) #41

#g1<-grep("HDL.",gwa$gwa)
#g2<-grep("LDL.",gwa$gwa)
#g3<-grep("VLDL.",gwa$gwa)
#g4<-grep("IDL.",gwa$gwa)

#g<-c(g1,g2,g3,g4) #112 traits
#gwa<-gwa[-g,]
#length(unique(gwa$gwa))
#[1] 140

#w<-which(gwa$gwas%in%c("Ace","AcAce","Ala","Alb","bOHBut","Bis.DB.ratio","Bis.FA.ratio","Cit","CH2.DB.ratio","CH2.in.FA","DB.in.FA","DHA","Est.C","Free.C","FAw3","FAw6","FAw67","FAw79S","Glc","Gln","Gly","Gp","His","Ile","LA","Lac","Leu","MUFA","Phe","Pyr","PC","otPUFA","Tyr","Tot.FA","Val","Serum.C","Serum.TG","SM","TotPG"))
#gwa<-gwa[-w,]
#length(unique(gwa$gwa)) #102

#w<-which(gwa$gwas%in%c("SDS"))
#gwa<-gwa[-w,]
#length(unique(gwa$gwa)) #101

#duplicated GWAS
#w<-which(gwa$gwas%in%c("CAD_2014","Alzheimers","ALS_lmm"))
#gwa<-gwa[-w,]

#length(unique(gwa$gwa)) #102

gwa %>%
    group_by(Annotation) %>%
    dplyr::summarize(Mean = mean(logOR, na.rm=TRUE))

# A tibble: 3 x 2
#  Annotation       Mean
#       <chr>      <dbl>
#1   cis only  0.9511548
#2  cis+trans -2.2809449
#3 trans only -8.8027932

# A tibble: 3 x 2
#Annotation        Mean
#<chr>       <dbl>
#  1   cis only   0.6871599
#2  cis+trans   0.7491646
#3 trans only -11.3796267

nthresh<-gwa %>%
    group_by(Phenotype) %>%
    dplyr::summarize(sum = sum(NThresh, na.rm=TRUE))

df<-data.frame(nthresh)
w<-which(df$sum==0)
length(w) #0

#w<-which(gwa$Phenotype%in%df[w,1]) #81
#gwa<-gwa[-w,] #traits
#length(unique(gwa$gwa)) #74

gwa %>%
    group_by(Annotation) %>%
    dplyr::summarize(Mean = mean(logOR, na.rm=TRUE))
# A tibble: 3 x 2
#Annotation        Mean
#<chr>       <dbl>
#  1   cis only   0.6871599
#2  cis+trans   0.7491646
#3 trans only -11.3796267
nthresh<-gwa %>%
    group_by(Phenotype) %>%
    dplyr::summarize(sum = sum(NThresh, na.rm=TRUE))

df<-data.frame(nthresh)
w<-which(df$sum<10)
length(w) #0
#w<-which(gwa$gwas%in%df[w,1])
#gwa<-gwa[-w,]#33
#length(unique(gwa$gwa)) #63

gwa$FDR<-p.adjust(gwa$Pvalue,method="BH")
write.table(gwa,"GARFIELD_gwa_traits.txt",sep="\t",row.names=F,col.names=T,quote=F)

gwa2<-gwa[which(gwa$FDR<0.05),]
table(gwa2$Annotation)

#cis only  cis+trans trans only
#14         19          1
gwa<-data.frame(gwa,gwa="mQTL SNPs amongst GWA SNPs")
p1<-ggplot(gwa, aes(x=logOR)) +
geom_density(aes(fill=Annotation), alpha=0.5, bw=0.1,adjust=3) +
scale_fill_brewer(type="qual") +
geom_vline(xintercept=0, linetype="dotted") +
facet_wrap(~gwa,nrow=1)+
theme(legend.position="bottom") +
theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
theme(strip.text.x = element_text(size = 16)) +
xlim(-6, 6) +
ylim(0, 1) +
labs(x="Enrichment (log OR)")

df2$annotation<-gsub("CpGs","DNAm sites",df2$annotation)
p2<-ggplot(df2, aes(x=logOR)) +
  geom_density(aes(fill=ciscat), alpha=0.5, bw=0.1) +
  scale_fill_brewer(type="qual") +
  labs(fill = "Annotation") +
  facet_wrap(~annotation,nrow=1)+
  geom_vline(xintercept=0, linetype="dotted") +
  theme(legend.position="bottom") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) +
  theme(strip.text.x = element_text(size = 16)) +
  labs(x="Enrichment (log OR)") +
  xlim(-3, 2.75) +
  ylim(0, 3)


pdf("Fig2.v2.pdf",height=10,width=16)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p2, vp = vplayout(1, 1:2))
print(p1, vp = vplayout(1, 3))
dev.off()

