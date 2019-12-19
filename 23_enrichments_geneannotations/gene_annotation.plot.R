#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)

r<-read.table(paste0(path,"mqtl_epigenetic/garfield.test.mqtl_epigenetic.out"),he=T)
r$cis_snp<-"All"
length(unique(r$Type))

cis<-read.table(paste0(path,"mqtl_cis_epigenetic/garfield.test.mqtl_cis_epigenetic.out"),he=T)
cis$cis_snp<-"cis only"
length(unique(cis$Type))

trans<-read.table(paste0(path,"mqtl_trans_epigenetic/garfield.test.mqtl_trans_epigenetic.out"),he=T)
trans$cis_snp<-"trans only"

amb<-read.table(paste0(path,"mqtl_ambivalent_epigenetic/garfield.test.mqtl_ambivalent_epigenetic.out"),he=T)
amb$cis_snp<-"cis+trans"

padj<-read.table(paste0(path,"mqtl_epigenetic/garfield.Meff.mqtl_epigenetic.out"),he=F)
padj<-padj[which(padj$V1=="Padj"),"V2"]

df<-rbind(r,cis,trans,amb)

df<-df[which(df$Category=="Genic"),]
df$logOddsRatio<-log(df$OR)

df2<-df[df$PThresh=="1e-14",]
nvar<-df2%>%group_by(cis_snp)%>%summarize(mean(NThresh))

nvar
# A tibble: 4 x 2
#  cis_snp    `mean(NThresh)`
#  <chr>                <dbl>
#1 All                 159270
#2 cis only            151498
#3 cis+trans            29637
#4 trans only            2641


#df2$cis_snp<-gsub("All",paste0("All (N=",nvar[1,2]," SNPs)"),df2$cis_snp)
#df2$cis_snp<-gsub("cis+trans",paste0("cis+trans (N=",nvar[2,2]," SNPs)"),df2$cis_snp)
#df2$cis_snp<-gsub("cis only",paste0("cis only (N=",nvar[3,2]," SNPs)"),df2$cis_snp)
#df2$cis_snp<-gsub("trans only",paste0("trans only (N=",nvar[4,2]," SNPs)"),df2$cis_snp)

df2$Type<-as.factor(df2$Type)
df2$Category<-as.factor(df2$Category)
cats<-levels(df2$Category)

#7 annotations for 4 categories
#pval_lim<-padj

pval_lim<-1e-3

w<-which(df2$Pvalue==0)
m<-min(df2$Pvalue[-w])
df2$Pvalue[w]<-m
length(unique(df2$Type)) #7

r.all2<-df2
r.all3<-r.all2[which(r.all2$cis_snp!="All"),]
names(r.all3)<-gsub("NAnnotThesh","NAnnotThresh",names(r.all3))
names(r.all3)<-gsub("Annotation","Gene_Annotation",names(r.all3))
names(r.all3)<-gsub("cis_snp","Annotation",names(r.all3))
names(r.all3)<-gsub("Type","Gene Annotation",names(r.all3))
r.all3<-r.all3[,c("Gene Annotation","OR","Pvalue","Beta","SE","CI95_lower","CI95_upper","NAnnotThresh","NAnnot","NThresh","N","Annotation")]

write.table(r.all3,"TableSXX_garfield_gene_annotation.txt",quote=F,row.names=F,col.names=T,sep="\t")
#State  OR  Pvalue  FDR_Pvalue  support b   c   d   CellType    Tissue  Annotation                                                                                                                                              

n<-which(names(df2)%in%c("Annotation","Celltype","Tissue","logOddsRatio"))
df3<-df2[which(df2$Category=="Genic"),-n]

n<-which(names(df3)%in%"cis_snp")
names(df3)[n]<-"Annotation"

p1<-ggplot(df3, aes(fill=Annotation, y=-log10(df3$Pvalue), x=Type)) + 
    geom_bar(position="dodge", stat="identity") +
    labs(x="",y="-log10 (Pvalue)") +
    theme(axis.text.x = element_blank())
    #scale_fill_brewer(type="qual")

m1<-min(round(df3$Beta,0))
m2<-max(round(df3$Beta,0))

p2<-ggplot(df3, aes(color=Annotation, y=Beta, x=Type)) +
 #scale_fill_brewer(type="qual") +
    #geom_bar(position="dodge", stat="identity") +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(x=Type,ymin=CI95_lower, ymax=CI95_upper), width=.2,position=position_dodge(.9)) +
    #scale_fill_brewer(type="qual") +
    labs(x="Gene annotation",y="log OR") +
    #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1:m),limits=c(0.2,m)) +
    scale_y_continuous(limits=c(-m2,m2)) +
    geom_hline(aes(yintercept = 0))
    #scale_fill_brewer(type="qual")
    #theme(axis.text.x=element_text(angle = 90, hjust = 10))    
#ggsave(p1,file="test.pdf",height=6,width=10)

pdf("gene_annotation_enrichment.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))

dev.off()
