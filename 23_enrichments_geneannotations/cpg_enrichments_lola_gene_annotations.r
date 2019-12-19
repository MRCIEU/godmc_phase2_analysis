#module add languages/R-3.5.1-ATLAS-gcc-6.1

library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(biovizBase)
library(meffil)
library(dplyr)
library(gridExtra)
library(grid)

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

setwd("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments")

##load regiondb
regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/",collection="7regions")

#
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
max(clumped[which(clumped$cis==TRUE),"pval"])
#1e-4
max(clumped[which(clumped$cis==FALSE),"pval"])
#5e-8
retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))

data=as.data.table(clumped)

#get 450k locations
#Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
Illumina450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Illumina450_dt=as.data.table(Illumina450[,1:3])
Illumina450_dt=as.data.table(Illumina450)
Illumina450_dt[,cpgID:=row.names(Illumina450),]
Illumina450_dt[,cpgstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
Illumina450_dt[,cpgend_pre:=ifelse(strand=="-",pos+500,pos+501),]

#collapse overlaps
gr_range = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre)))
gr_cpg = with(Illumina450_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

Illumina450_dt[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,cpgend:=end(gr_range[overlap_red$subjectHit])]
Illumina450_dt[,NsubjectHits:=overlap_red$NsubjectHits]

create_grs = function(data,selector){
  fg_cpg_list=GRangesList()
  for (sel in selector){
    fg_cpg=with(unique(data[eval(parse(text=sel)),c("cpgchr","cpgstart","cpgend"),with=FALSE]),GRanges(seqnames=cpgchr,ranges=IRanges(cpgstart,cpgend)))
    fg_cpg_list[[paste0(sel,"_cpg ",length(fg_cpg)," regions")]]=fg_cpg
    
  }
 
  bg_list=GRangesList()
  bg_list[["cpg"]]=with(unique(data[,c("cpgchr","cpgstart","cpgend"),with=FALSE]),GRanges(seqnames=cpgchr,ranges=IRanges(cpgstart,cpgend)))
    
  return(list(fg_cpg=fg_cpg_list,bg=bg_list))
}

Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==cpgpos,])

data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]

data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'"))

test<-unique(data.frame(data$cpg,data$cpg_cis))
table(test[,2])
#ambivalent      FALSE       TRUE 
#     28023       4267     157812 


hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)
#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~3000)
Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isGoDMC:=ifelse(cpgID%in%data$cpg,TRUE,FALSE),]
Illumina450_dt<-Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]

df<-unique(data.frame(cpg=data$cpg,cpg_cis=data$cpg_cis))
m<-match(Illumina450_dt$cpgID,data$cpg)
df<-data[m,"cpg_cis"]
Illumina450_dt$cpg_cis<-df

Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$quantiles<-paste(Illumina450_dt$CpG_freqquantile,Illumina450_dt$GC_freqquantile)

#ALL
t<-table(Illumina450_dt$quantiles,Illumina450_dt$isGoDMC)
w<-which(t[,2]<5)
pr<-data.frame(t[-w,1]/t[-w,2])
m<-min(pr)*t[,2]

controllist<-list()
for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[i,1]>5&t[i,2]>0){
subgroup<-Illumina450_dt[which(Illumina450_dt$quantiles==row.names(t)[i]&Illumina450_dt$isGoDMC==FALSE& Illumina450_dt$chr!="chrY"),]
id<-subgroup[sample(nrow(subgroup), size=round(m[i],0), replace=FALSE),"cpgID"]
controllist[[i]]<-id
}}

bg.matched<-do.call("rbind",controllist)
Illumina450_bg_matched<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched$cpgID),]
Illumina450_godmc<-Illumina450_dt[which(Illumina450_dt$cpgID%in%data$cpg),]


bg_matched<-rbind(Illumina450_bg_matched,Illumina450_godmc)

GoDMC_cpg_gr=unique(with(data,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

#use all for background
Illumina450_bg=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

lola_res0=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB, cores=5)
lola_res0$logOddsRatio<-log(lola_res0$oddsRatio)
#plotLOLA(locResults_all=lola_res0,plot_pref="cpg_gene450kbg",height=10,width=18)

#Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
#lola_res0_matched=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)
Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

cpg_bg_gr_matched=unique(c(Illumina450_bg_matched,GoDMC_cpg_gr))

lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
plotLOLA(locResults_all=lola_res0_matched,plot_pref="cpg_corebg_matched",height=10,width=18)
save(lola_res0_matched,lola_res0,file="../results/enrichments/lola_gene_mqtlcpg.rdata")

##ambivalent/FALSE/TRUE

w<-which(is.na(Illumina450_dt$cpg_cis))
Illumina450_dt$cpg_cis[w]<-"0bg"
n<-names(table(Illumina450_dt$cpg_cis))[-1]

t<-table(Illumina450_dt$quantiles,Illumina450_dt$cpg_cis)
bg.matched.subset<-list()
Illumina450_bg_matchedcis<-list()
res_all<-list()
l<-list()

data.frame(t[,2]/sum(t[,2]),t[,3]/sum(t[,3]),t[,4]/sum(t[,4]))
a<-apply(t,2,sum)
#identify column with most annotations
col<-apply(t[,-1],2,sum)
col<-which(col==max(col))+1

for (j in 1:length(n)){
cat(n[j],"\n")
#w<-which(t[,col]<5)
#pr<-data.frame(t[-w,1]/t[-w,col])
w<-which(t[,(j+1)]<5)
pr<-data.frame(t[-w,1]/t[-w,(j+1)])
m<-min(pr)*t[,(j+1)]

controllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[i,1]>5&t[i,(j+1)]>0){
subgroup<-Illumina450_dt[which(Illumina450_dt$quantiles==row.names(t)[i]&Illumina450_dt$isGoDMC==FALSE),]
id<-subgroup[sample(nrow(subgroup), size=round(m[i],0), replace=FALSE),"cpgID"]
controllist[[i]]<-id
}}

bg.matched.subset[[j]]<-do.call("rbind",controllist)
Illumina450_bg_matchedcis[[j]]<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched.subset[[j]]$cpgID),]
Illumina450_bg_matchedcis[[j]]<-unique(with(Illumina450_bg_matchedcis[[j]],GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
data_subset<-data[which(data$cpg_cis==n[j]),]
GoDMC_cpg_gr=unique(with(data_subset,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
l[[j]]<-length(GoDMC_cpg_gr)


cpg_bg_gr_matched=unique(c(Illumina450_bg_matchedcis[[j]],GoDMC_cpg_gr))
lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
lola_res0_matched$userSet<-j
res_all<-rbind(res_all,lola_res0_matched)

}
res_all$Annotation<-res_all$userSet
res_all$Annotation<-gsub("ambivalent (10941 regions)","cis+trans",res_all$Annotation,fixed=T)
res_all$Annotation<-gsub("cis_only (107004 regions)","cis only",res_all$Annotation,fixed=T)
res_all$Annotation<-gsub("trans_only (3864 regions)","trans only",res_all$Annotation,fixed=T)

save(res_all,file="../results/enrichments/lola_gene_mqtlcpg_cis.rdata")

res_all[,Pvalue:=10^(-pValueLog),]
 
res_all[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
  res_all[,mlog10p.adjust:=-log10(p.adjust),]
  res_all$description<-gsub("hg19_genes_","",res_all$description)
  res_all$description<-gsub("hg19_","",res_all$description)


#res_all2<-res_all[which(res_all$p.adjust<0.001),]
#dim(res_all2)
res_all2<-res_all
res_all3<-res_all2[,c("description","oddsRatio","Pvalue","p.adjust","support","b","c","d","Annotation")]
names(res_all3)<-c("Gene Annotation","OR","Pvalue","FDR_Pvalue","support","b","c","d","Annotation")                                                   

write.table(res_all3,"TableSXX_LOLA_gene_annotation_updated.txt",sep="\t",quote=F,row.names=F,col.names=T)


plotLOLA_OR=function(locResults_all,plot_pref,height=35,width=18){
  
  ##process LOLA results
  #pval_lim=0.001
  
  #locResults=process_LOLA(LOLA_res=locResults_all,cellType_conversions=cellType_conversions)
  
  #sub=locResults
  #sub[,signif:=any(p.adjust<=pval_lim),by="seg_explanation"]
  #sub=sub[signif==TRUE]
  #sub[,tissue_count:=paste0(tissue," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",tissue_count_all),by=c("tissue")]
  
   #changed form exp() to 10^ to accomodate change in LOLA
  locResults<-locResults_all
   
   ##bubble plots
  m<-max(locResults$oddsRatio)+1
  #if(locResults$description[1]=="25states"){
  pdf(paste0(plot_pref,"_All_OR.pdf"),height=height,width=width+2)
  pl3=ggplot(locResults,aes(x=description,y=oddsRatio))+
  geom_hline(yintercept=1, linetype="dotted")+
  geom_point(position=position_dodge(.9)) +
  geom_errorbar(aes(x=description,ymin=conf_down, ymax=conf_up), width=.2,position=position_dodge(.9)) +
   
 facet_wrap(~Annotation,scale="free_x",ncol=1)+
  #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  #scale_size(range=c(1,8))+
  #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  #ylim(0.1,4)+
  #guides(colour = guide_legend(override.aes = list(size=10,ncol=10))) +
  #guides(fill = guide_legend(ncol=20))+
  #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,6),limits=c(0.2,6)) +
  #theme(legend.text=element_text(size=12)) +
  labs(y="Odds ratio (log scale)",x="Gene annotation")
  #scale_fill_brewer(type="qual")
  
  print(pl3)
  dev.off()
} 
#}

plotLOLA_OR(locResults_all=res_all,plot_pref="cpg_gene_matched_cis_OR",height=10,width=18)

res_all$description <- factor(x = res_all$description, levels = c("cpg_shelves","cpg_shores","cpg_islands","cpg_inter","1to5kb","promoters","5UTRs","exons","introns","intronexonboundaries","3UTRs","intergenic"))
                                       
p1<-ggplot(res_all, aes(fill=Annotation, y=-log10(res_all$pvalue), x=description)) + 
    geom_bar(position="dodge", stat="identity") +
    labs(x="",y="-log10 (Pvalue)") +
    theme(axis.text.x = element_blank())
    #scale_fill_brewer(type="qual")

m1<-min(round(res_all$logOddsRatio,0))
m2<-max(round(res_all$logOddsRatio,0))

p2<-ggplot(res_all, aes(color=Annotation, y=logOddsRatio, x=description)) +
 #scale_fill_brewer(type="qual") +
    #geom_bar(position="dodge", stat="identity") +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(x=description,ymin=log(res_all$conf_down), ymax=log(res_all$conf_up)), width=.2,position=position_dodge(.9)) +
    #scale_fill_brewer(type="qual") +
    labs(x="Gene annotation",y="log OR") +
    #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1:m),limits=c(0.2,m)) +
    scale_y_continuous(limits=c(-m2,m2)) +
    geom_hline(aes(yintercept = 0))
    #scale_fill_brewer(type="qual")
    #theme(axis.text.x=element_text(angle = 90, hjust = 10))    
#ggsave(p1,file="test.pdf",height=6,width=10)

pdf("gene_annotation_cpg_enrichment.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))

dev.off()

