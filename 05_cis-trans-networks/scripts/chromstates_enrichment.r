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
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

##Functions
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

process_LOLA = function (LOLA_res, collections=c("25states"),cellType_conversions){
  
  LOLA_res=LOLA_res[!is.na(userSet)]
  
  LOLA_res=LOLA_res[collection%in%collections]

  #changed form exp() to 10^ to accomodate change in LOLA   
  w<-which(LOLA_res$pValueLog=="Inf")
  if(length(w)>0){
  LOLA_res$pValueLog[w]<--log10(2.5e-320) }
  LOLA_res[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
  LOLA_res[,mlog10p.adjust:=-log10(p.adjust),]
  LOLA_res[,tissue_count_all:=length(unique(filename[!is.na(filename)])),by=c("tissue")]
  LOLA_res[,target:=toupper(sub("-","",unlist(lapply(antibody,function(x){spl=unlist(strsplit(x,"_|eGFP-"));spl[spl!=""][1]})))),]
  
  return(LOLA_res)  
  
}

plotLOLA_OR=function(locResults_all,plot_pref,height=350,width=18){
  
  ##process LOLA results
  pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=locResults_all,cellType_conversions=cellType_conversions)
  
  sub=locResults
  sub[,signif:=any(p.adjust<=pval_lim),by="seg_explanation"]
  sub=sub[signif==TRUE]
  sub=sub[p.adjust!=1]
  #sub[,tissue_count:=paste0(tissue," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",tissue_count_all),by=c("tissue")]
  
   ##bubble plots
  m<-max(locResults$oddsRatio)+1
  if(locResults$description[1]=="25states"){
  pdf(paste0(plot_pref,"_All_OR.pdf"),height=height,width=width+2)
  pl3=ggplot(locResults,aes(x=seg_explanation,y=oddsRatio))+
  geom_hline(yintercept=1, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=tissue,size=pValueLog))+
  facet_wrap(~cluster,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  #ylim(0.1,4)+
  #guides(colour = guide_legend(override.aes = list(size=10,ncol=10))) +
  guides(fill = guide_legend(ncol=20))+
  scale_y_continuous(trans = 'log10',breaks=c(.1,.2,.3,.4,.5,1,2,3,4,5,10,20,50),limits=c(0.1,50)) +
  theme(legend.text=element_text(size=12)) +
  labs(y="Odds ratio (log scale)",x="State") +
  scale_fill_brewer(type="qual")
  
  print(pl3)
  dev.off()
} 
}

core_communities_cpg_tophits$userSet<-factor(core_communities_cpg_tophits$userSet,levels=clusters)
  p1<-ggplot(core_communities_cpg_tophits,aes(x=factor(userSet),y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=antibody2))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
      guides(fill = guide_legend(nrow=2))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Transcription factor",fill="TFBS")
    ggsave(p1,file=paste0("../images/epigenetic_TFBSbycommunity_OR.pdf"),height=10,width=16)

##load regiondb
regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states",collection="25states")


##load cell type conversion and colors
cellType_conversions=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/CellTypes.tsv",drop="collection")
colors=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/color.tsv")

#setwd("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments")
load("../results/graph.rdata")
dim(mem)
#5109 cpgs 
length(unique(mem$cluster))
#[1] 1615
sum(table(mem$cluster) >= 10)
#56
cl<-data.frame(table(mem$cluster))
cl<-cl[which(cl$Freq>=10),1]
mem2<-mem[which(mem$cluster%in%cl),]
length(table(mem2$cluster))
#56
y<-meffil.get.features("450k")
m<-match(mem2$cpg,y$name)
mem2<-data.frame(mem2,cpgchr=y[m,"chromosome"],cpgpos=y[m,"position"])
data=as.data.table(mem2)

#get 450k locations
Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
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

##prepare data: merge CpGs that are in proximity to eachother to avoid infalting the results, 1kb around cp
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==cpgpos,])

data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]
#data[,cpg_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]

####run with internal background
#grs_mqtl_effect=create_grs(data=data,selector=c("cpg_change=='mqtl_effect>0'","cpg_change=='mqtl_effect<0'","cpg_change=='ambivalent'","snp_change=='mqtl_effect>0'","snp_change=='mqtl_effect<0'","snp_change=='ambivalent'"))
#produceLOLA_plots_intern(grs=grs_mqtl_effect,plot_pref="mqtl_effect",height=35,width=18)
#data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
#grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'"))

data[,cluster,by=c("cpgchr","cpgstart","cpgend")]
grs_cis=create_grs(data=data,selector=unique(data$cluster))
length(grs_cis$fg_cpg) #56
length(grs_cis$bg) #1

hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)
#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~3000)
Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isGoDMC:=ifelse(cpgID%in%data$cpg,TRUE,FALSE),]
Illumina450_dt<-Illumina450_dt[Illumina450_dt$cpgID%in%unique(mem$cpg),]

table(Illumina450_dt$isGoDMC)

#FALSE  TRUE 
# 3779  1330 

#df<-unique(data.frame(cpg=data$cpg,cpg_cis=data$cpg_cis))
m<-match(Illumina450_dt$cpgID,data$cpg)
df<-data[m,"cluster"]
Illumina450_dt$cluster<-df

#Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
#Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)

Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)
Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)

Illumina450_dt$quantiles<-paste(Illumina450_dt$CpG_freqquantile,Illumina450_dt$GC_freqquantile)

##

w<-which(is.na(Illumina450_dt$cluster))
Illumina450_dt$cluster[w]<-"0bg"
n<-names(table(Illumina450_dt$cluster))[-1]

t<-table(Illumina450_dt$quantiles,Illumina450_dt$cluster)
t2<-table(Illumina450_dt$quantiles,Illumina450_dt$isGoDMC)

bg.matched.subset<-list()
Illumina450_bg_matchedcis<-list()
res_all<-list()
l<-list()

#number of clusters
for (j in 1:length(n)){
cat(n[j],"\n") #cluster id
w<-which(t[,13]<5) #clusters with less than 5 cpgs
pr<-data.frame(t[-w,1]/t[-w,13])
m<-min(pr)*t[,(j+1)]

controllist<-list()

#nrow(t) is number of bins
for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[i,1]>5&t[i,2]>0){
subgroup<-Illumina450_dt[which(Illumina450_dt$quantiles==row.names(t)[i]&Illumina450_dt$cluster!=n[j]),]
id<-subgroup[sample(nrow(subgroup), size=round(m[i],0), replace=FALSE),"cpgID"]
controllist[[i]]<-id
}}

bg.matched.subset[[j]]<-do.call("rbind",controllist)
Illumina450_bg_matchedcis[[j]]<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched.subset[[j]]$cpgID),]
Illumina450_bg_matchedcis[[j]]<-unique(with(Illumina450_bg_matchedcis[[j]],GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
data_subset<-data[which(data$cluster==n[j]),]
GoDMC_cpg_gr=unique(with(data_subset,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
l[[j]]<-length(GoDMC_cpg_gr)


cpg_bg_gr_matched=unique(c(Illumina450_bg_matchedcis[[j]],GoDMC_cpg_gr))
lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
lola_res0_matched$userSet<-j
res_all<-rbind(res_all,lola_res0_matched)

}
cl_id<-data.frame(userSet=1:length(n),cluster=colnames(t)[-1])
m<-match(res_all$userSet,cl_id$userSet)
res_all$cluster<-cl_id[m,"cluster"]

spl<-do.call("rbind",strsplit(res_all$filename,split="_"))
res_all$dataSource<-spl[,1]
res_all$seg_code<-gsub(".bed","",spl[,5])

state<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T)
state$STATE<-paste0("E",state$STATE)
m<-match(res_all$seg_code,state$STATE)
res_all$seg_explanation<-state[m,"NO."]
res_all$seg_explanation<-gsub("prime","'",res_all$seg_explanation)
#res_all$seg_explanation<-as.character(res_all$seg_explanation)
w<-which(is.na(m))
res_all<-res_all[-w,]
length(unique(res_all$seg_code))
#25
tiss<-read.table("../../07_enrichments/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
m<-match(res_all$dataSource,as.character(tiss$ID))
res_all$tissue<-tiss[m,"Tissue"]
res_all$tissue<-gsub("GI_","",res_all$tissue)
res_all$tissue<-tolower(res_all$tissue)


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


res_all$tissue<-sapply(res_all$tissue, simpleCap)
res_all$tissue<-gsub("Esc","ESC",res_all$tissue)
res_all$tissue<-gsub("Ipsc","iPSC",res_all$tissue)
res_all$Annotation<-res_all$userSet
res_all$Annotation<-gsub("ambivalent (10941 regions)","cis+trans",res_all$Annotation,fixed=T)
res_all$Annotation<-gsub("cis_only (107004 regions)","cis only",res_all$Annotation,fixed=T)
res_all$Annotation<-gsub("trans_only (3864 regions)","trans only",res_all$Annotation,fixed=T)
tiss<-read.table("~/repo/godmc_phase2_analysis/07_enrichments/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
m<-match(res_all$dataSource,tiss$ID)
res_all$cellType<-tiss$Name[m]

w<-which(res_all$cluster=="2")
res_all$Pval<-10^(-res_all$pValueLog)
plotLOLA_OR(locResults_all=res_all,plot_pref="cpg_extbg_chromstates_matched_cis_OR",height=350,width=18)

w<-which(res_all$cluster%in%c("1","2","6","10","14"))
res_all_gwa<-res_all[w,]
max(res_all_gwa$oddsRatio[which(res_all_gwa$oddsRatio!="Inf")])
#[1] 23.88923
plotLOLA_OR(locResults_all=data.table(res_all_gwa),plot_pref="cpg_extbg_chromstates_matched_cis_OR_gwacluster",height=30,width=18)

res_all<-res_all[which(res_all$logOddsRatio!="-Inf"),]
res_all[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
res_all[,Pvalue:=10^(-pValueLog)]

res_all2<-res_all[which(res_all$p.adjust<0.001),]

res_all3<-res_all2[,c("cluster","seg_explanation","oddsRatio","Pvalue","p.adjust","support","b","c","d","cellType","tissue")]
names(res_all3)<-c("Community","State","OR","Pvalue","FDR_Pvalue","support","b","c","d","CellType","Tissue")                                                   

write.table(res_all3,"TableSXX_LOLA_segmentationstate_updated.txt",sep="\t",quote=F,row.names=F,col.names=T)


dim(res_all2)
unique(res_all2$cluster)
#[1] 2  76
unique(res_all2[,c("cluster","seg_explanation")])
#   cluster seg_explanation
#1:       2           EnhA1
#2:      76          ReprPC
#3:      76         PromBiv

w<-which(res_all$cluster%in%c("2"))
plotLOLA_OR(locResults_all=res_all[w,],plot_pref="cpg_extbg_chromstates_matched_cis_OR_2",height=10,width=18)
w<-which(res_all$cluster%in%c("76"))
plotLOLA_OR(locResults_all=res_all[w,],plot_pref="cpg_extbg_chromstates_matched_cis_OR_76",height=10,width=18)
save(res_all,file="../results/lola_chromstates_communities.rdata")
