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

#r<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19/roadmap_epigenomics/index.tmp",he=T)
#table(r$dataType)

#  dnase histone 
#    131     979 


##load cell type conversion and colors
cellType_conversions=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/CellTypes.tsv",drop="collection")
colors=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/color.tsv")

setwd("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments")

##load regiondb
regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states",collection="25states")

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

plotLOLA=function(locResults_all,plot_pref,height=35,width=18){
  
  ##process LOLA results
  pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=locResults_all,cellType_conversions=cellType_conversions)
  
  sub=locResults
  sub[,signif:=any(p.adjust<=pval_lim),by="seg_explanation"]
  sub=sub[signif==TRUE]
  sub[,tissue_count:=paste0(tissue," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",tissue_count_all),by=c("tissue")]
  
  ##summary plot
  sub_sum=sub[,list(mean_mlog10p.adjust=mean(mlog10p.adjust[p.adjust<=pval_lim&mlog10p.adjust!=Inf]),count=length(userSet[p.adjust<=pval_lim])),by="userSet"]
  
  pdf(paste0(plot_pref,"_summary.pdf"),height=height/7,width=10)
  pl=ggplot(sub_sum,aes(x=userSet,y=count))+geom_bar(stat="identity",aes(fill=mean_mlog10p.adjust))+ylab("Enriched region sets")+coord_flip()+scale_fill_gradient(low="blue",high="red")
  print(pl)
  dev.off()
  
  ##bubble plots
  if(sub$description[1]=="25states"){
  pdf(paste0(plot_pref,"_All.pdf"),height=height,width=width+2)
  pl3=ggplot(sub,aes(x=seg_explanation,y=-log10(p.adjust),size=logOddsRatio,fill=tissue))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  ylim(0,350)+
  guides(fill = guide_legend(ncol=10))+
  theme(legend.text=element_text(size=10))
  
  print(pl3)
  dev.off()
} 
}

plotLOLA_OR=function(locResults_all,plot_pref,height=35,width=18){
  
  ##process LOLA results
  #pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=locResults_all,cellType_conversions=cellType_conversions)
  
  #sub=locResults
  #sub[,signif:=any(p.adjust<=pval_lim),by="seg_explanation"]
  #sub=sub[signif==TRUE]
  #sub[,tissue_count:=paste0(tissue," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",tissue_count_all),by=c("tissue")]
  
   ##bubble plots
  m<-max(locResults$oddsRatio)+1
  if(locResults$description[1]=="25states"){
  pdf(paste0(plot_pref,"_All_OR.pdf"),height=height,width=width+2)
  pl3=ggplot(locResults,aes(x=seg_explanation,y=oddsRatio))+
  geom_hline(yintercept=1, linetype="dotted")+
  geom_jitter(width = 0.2, aes(colour=tissue,size=pValueLog))+
  facet_wrap(~Annotation,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,8))+
  #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+
  #ylim(0.1,4)+
  #guides(colour = guide_legend(override.aes = list(size=10,ncol=10))) +
  guides(fill = guide_legend(ncol=20))+
  scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,6),limits=c(0.2,6)) +
  theme(legend.text=element_text(size=12)) +
  labs(y="Odds ratio (log scale)",x="State") +
  scale_fill_brewer(type="qual")
  
  print(pl3)
  dev.off()
} 
}

#produceLOLA_plots_intern = function(grs,plot_pref,height=35,width=18,recreate=FALSE){

  #run LOLA cpgs
 # simpleCache(cacheName=paste0(plot_pref,"_cpgs"),instruction="runLOLA(fg_cpg, bg_cpg, regionDB0, cores=5)",cacheDir=getwd(),recreate=recreate,assignToVariable="locResults_cpgs",buildEnvir=c(fg_cpg=grs$fg_cpg,bg_cpg=grs$bg$cpg))
   
  #combine and process
 # locResults_all=locResults_cpgs
 # locResults_all$logOddsRatio<-log(locResults_all$oddsRatio)
 # plotLOLA(locResults_all=locResults_all,plot_pref=plot_pref,height=height,width=width)
#}

##prepare GoDMC data: merge CpGs and SNPs that are in proximity to eachother to avoid infalting the results, 1kb around cp
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==cpgpos,])

data[,cpgchr:=gsub("23","X",cpgchr),]
data[,cpgchr:=gsub("24","Y",cpgchr),]
#data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpg")]
#test<-unique(data.frame(data$cpg,data$cpg_cis))
#table(test[,2])

#ambivalent      FALSE       TRUE 
#     11902       7214     170986 
data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
test<-unique(data.frame(data$cpg,data$cpg_cis))
table(test[,2])
#ambivalent      FALSE       TRUE 
#     28023       4267     157812

grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'"))

####run with external background for CpGs

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

#plot CG and CpG frequency for GoDMC cpgs and background 
pdf("./images/compare_seq_properties.pdf",height=3,width=4)
ggplot(Illumina450_dt,aes(x=GC_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()
dev.off()

pdf("./images/compare_seq_properties_matched.pdf",height=3,width=4)
ggplot(bg_matched,aes(x=GC_freq,col=isGoDMC))+geom_density()
ggplot(bg_matched,aes(x=GC_freq))+geom_density()
ggplot(bg_matched,aes(x=CpG_freq,col=isGoDMC))+geom_density()
ggplot(bg_matched,aes(x=CpG_freq))+geom_density()
dev.off()


#use all for background
Illumina450_bg=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

#lola_res0=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB, cores=5)
#lola_res0$logOddsRatio<-log(lola_res0$oddsRatio)



#plotLOLA(locResults_all=lola_res0,plot_pref="cpg_ext450kbg",height=10,width=18)

#Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
#lola_res0_matched=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)
Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

cpg_bg_gr_matched=unique(c(Illumina450_bg_matched,GoDMC_cpg_gr))

lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
save(lola_res0_matched,file="../results/enrichments/lola_chromstates_mqtlcpg.rdata")

#subset by dnase and histonemark
#w<-which(is.na(lola_res0_matched$antibody)&lola_res0_matched$collection=="roadmap_epigenomics")
#lola_res0_matched[w,]

spl<-do.call("rbind",strsplit(lola_res0_matched$filename,split="_"))
lola_res0_matched$dataSource<-spl[,1]
lola_res0_matched$seg_code<-gsub(".bed","",spl[,5])

state<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T)
state$STATE<-paste0("E",state$STATE)
m<-match(lola_res0_matched$seg_code,state$STATE)
lola_res0_matched$seg_explanation<-state[m,"NO."]
#lola_res0_matched$seg_explanation<-as.character(lola_res0_matched$seg_explanation)
#lola_res0_matched$seg_explanation[which(is.na(m))]<-"expanded"

w<-which(is.na(m))
lola_res0_matched<-lola_res0_matched[-w,]

tiss<-read.table("jul2013.roadmapData_tissues.txt",sep="\t",he=T)
m<-match(lola_res0_matched$dataSource,as.character(tiss$ID))
lola_res0_matched$tissue<-tiss[m,"Tissue"]

plotLOLA(locResults_all=lola_res0_matched,plot_pref="cpg_extbg_chromstates_matched",height=10,width=18)
save(lola_res0_matched,file="../results/enrichments/lola_chromstates_mqtlcpg.rdata")

#save(lola_res0_matched,lola_res0,file="../results/lola_ext_mqtlcpg.rdata")

##ambivalent/FALSE/TRUE

w<-which(is.na(Illumina450_dt$cpg_cis))
Illumina450_dt$cpg_cis[w]<-"0bg"
n<-names(table(Illumina450_dt$cpg_cis))[-1]
#[1] "ambivalent" "FALSE"      "TRUE"      

t<-table(Illumina450_dt$quantiles,Illumina450_dt$cpg_cis)
bg.matched.subset<-list()
Illumina450_bg_matchedcis<-list()
res_all<-list()
l<-list()
#proportion for each bin are different for each annotation category
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

w<-which(res_all$userSet==1)
res_all$userSet[w]<-paste0(n[1]," (", l[1]," regions)")
w<-which(res_all$userSet==2)
res_all$userSet[w]<-paste0("trans_only (", l[2]," regions)")
w<-which(res_all$userSet==3)
res_all$userSet[w]<-paste0("cis_only (", l[3]," regions)")

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
#res_all$seg_explanation[which(is.na(m))]<-"expanded"


plotLOLA(locResults_all=res_all,plot_pref="cpg_extbg_chromstates_matched_cis_updated",height=10,width=18)


save(res_all,file="../results/enrichments/lola_chromstates_mqtlcpg_cis_updated.rdata")

tiss<-read.table("jul2013.roadmapData_tissues.txt",sep="\t",he=T)
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



plotLOLA_OR(locResults_all=res_all,plot_pref="cpg_extbg_chromstates_matched_cis_OR_updated",height=10,width=18)

res_all<-res_all[which(res_all$logOddsRatio!="-Inf"),]
res_all[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
res_all[,Pvalue:=10^(-pValueLog)]

res_all2<-res_all[which(res_all$p.adjust<0.001),]
dim(res_all2)

res_all3<-res_all2[,c("seg_explanation","oddsRatio","Pvalue","p.adjust","support","b","c","d","cellType","tissue","Annotation")]
names(res_all3)<-c("State","OR","Pvalue","FDR_Pvalue","support","b","c","d","CellType","Tissue","Annotation")                                                   

write.table(res_all3,"TableSXX_LOLA_segmentationstate_updated.txt",sep="\t",quote=F,row.names=F,col.names=T)

ORs<-res_all %>%
  group_by(seg_explanation,userSet) %>%
  summarise(mean = mean(oddsRatio),meanPadj=mean(p.adjust))

ORs_cis<-ORs %>%
  group_by(userSet) %>%
  summarise(mean = mean(mean),meanPadj=mean(meanPadj))

####
notblood<-res_all[which(res_all$tissue!="Blood"),]
notblood<-notblood[order(notblood$p.adjust,decreasing=F),]

ORs_notblood<-notblood %>%
  group_by(seg_explanation,userSet) %>%
  summarise(mean = mean(oddsRatio),meanPadj=mean(p.adjust))
names(ORs_notblood)<-paste0("notblood.",names(ORs_notblood))

####
blood<-res_all[which(res_all$tissue=="Blood"),]
blood<-blood[order(blood$p.adjust,decreasing=F),]

ORs_blood<-blood %>%
  group_by(seg_explanation,userSet) %>%
  summarise(mean = mean(oddsRatio),meanPadj=mean(p.adjust))
names(ORs_blood)<-paste0("blood.",names(ORs_blood))
####

ORs_all<-data.frame(ORs,ORs_blood[,3:4],ORs_notblood[,3:4])
ORs_all<-data.frame(ORs_all,OR_bloodminusother=(ORs_all$blood.mean-ORs_all$notblood.mean), ORfold=(ORs_all$blood.mean/ORs_all$notblood.mean)) #positive is stronger in blood
mean(ORs_all[which(ORs_all$OR_bloodminusother>0&ORs_all$meanPadj<0.001),"ORfold"]) #[1] 1.122
ORs_all[which(ORs_all$OR_bloodminusother>0&ORs_all$meanPadj<0.001&ORs_all$userSet=="cis_only (107004 regions)"),]
ORs_all<-ORs_all[which(ORs_all$userSet=="cis_only (107004 regions)"),]

w<-which(res_all$tissue=="Blood")
res_all$blood<-"0-Not in Blood"
res_all$blood[w]<-"Blood"

w1<-which(res_all$Annotation=="cis only")
w2<-which(res_all$Annotation=="cis+trans")
w3<-which(res_all$Annotation=="trans only")
summary(lm(res_all$oddsRatio~res_all$blood))
summary(lm(res_all$oddsRatio[w1]~res_all$blood[w1]))
summary(lm(res_all$oddsRatio[w2]~res_all$blood[w2]))
summary(lm(res_all$oddsRatio[w3]~res_all$blood[w3]))

enh<-res_all[grep("Enh",res_all$seg_explanation),]
w1<-which(enh$Annotation=="cis only")
w2<-which(enh$Annotation=="cis+trans")
w3<-which(enh$Annotation=="trans only")
summary(lm(enh$oddsRatio~enh$blood))
summary(lm(enh$oddsRatio[w1]~enh$blood[w1]))
summary(lm(enh$oddsRatio[w2]~enh$blood[w2]))
summary(lm(enh$oddsRatio[w3]~enh$blood[w3]))

prom<-res_all[grep("Prom",res_all$seg_explanation),]
w1<-which(prom$Annotation=="cis only")
w2<-which(prom$Annotation=="cis+trans")
w3<-which(prom$Annotation=="trans only")
summary(lm(prom$oddsRatio~prom$blood))
summary(lm(prom$oddsRatio[w1]~prom$blood[w1]))
summary(lm(prom$oddsRatio[w2]~prom$blood[w2]))
summary(lm(prom$oddsRatio[w3]~prom$blood[w3]))

dnase<-res_all[grep("DNase",res_all$seg_explanation),]
summary(lm(dnase$oddsRatio~dnase$blood))
tx<-res_all[grep("Tx",res_all$seg_explanation),]
summary(lm(tx$oddsRatio~tx$blood))

#cis only  cis+trans trans only 
#      3175       3175       3175 


ORs<-data.frame(ORs)
o<-order(ORs$meanPadj)
ORs<-ORs[o,]

ORs1<-ORs[which(ORs$seg_explanation%in%c("PromD1","PromD2","PromBiv","PromU","PromP")),]
ORs1<-ORs1[which(ORs1$userSet=="cis_only (107004 regions)"),]
mean(ORs1$mean)
#1.313537

ORs1<-ORs[which(ORs$seg_explanation%in%c("PromD1","PromD2","PromBiv","PromU","PromP")),]
ORs1<-ORs1[which(ORs1$userSet=="ambivalent (10941 regions)"),]
mean(ORs1$mean)
#1.407478

ORs1<-ORs[which(ORs$seg_explanation%in%c("PromD1","PromD2","PromBiv","PromU","PromP")),]
ORs1<-ORs1[which(ORs1$userSet=="trans_only (3864 regions)"),]
mean(ORs1$mean)
#[1] 1.389386

ORs1<-ORs[which(ORs$seg_explanation%in%c("EnhW1","EnhW2","EnhAF","EnhA2","EnhA1","EnhAc")),]
ORs1<-ORs1[which(ORs1$userSet=="cis_only (107004 regions)"),]
mean(ORs1$mean)
#[1] 1.371817

ORs1<-ORs[which(ORs$seg_explanation%in%c("EnhW1","EnhW2","EnhAF","EnhA2","EnhA1","EnhAc")),]
ORs1<-ORs1[which(ORs1$userSet=="ambivalent (10941 regions)"),]
mean(ORs1$mean)
#1.648661

ORs1<-ORs[which(ORs$seg_explanation%in%c("EnhW1","EnhW2","EnhAF","EnhA2","EnhA1","EnhAc")),]
ORs1<-ORs1[which(ORs1$userSet=="trans_only (3864 regions)"),]
mean(ORs1$mean)
#[1] 0.9716326

ORs1<-ORs[which(ORs$seg_explanation%in%c("TssA")),]

write.table(res_all,"../results/enrichments/lola_chromstates_mqtlcpg_cis_updated.txt",sep="\t",col.names=T,quote=F,row.names=F)
res_all<-res_all[which(res_all$p.adjust<0.001),]

#
trans<-res_all[res_all$userSet=="trans_only (3864 regions)",]
df<-unique(data.frame(amb$seg_explanation,amb$tissue))
df<-data.frame(table(df[,1]))
dim(df) #127
dim(df[df$Freq>3,]) #12

amb<-unique(res_all[userSet=="ambivalent (10941 regions)",])
df<-unique(data.frame(amb$seg_explanation,amb$tissue))
df<-data.frame(table(df[,1]))
dim(df) #127
dim(df[df$Freq>3,]) #12

cis<-unique(res_all[userSet=="cis_only (107004 regions)",])
df<-unique(data.frame(cis$seg_explanation,cis$tissue))
df<-data.frame(table(df[,1]))
dim(df) #127
dim(df[df$Freq>3,]) #12

sub_blood<-sub[sub$tissue=="BLOOD",]

amb<-sub_blood[sub_blood$userSet=="ambivalent (10941 regions)",]
amb<-amb[amb$p.adjust<0.05&amb$oddsRatio>0,]
table(amb$seg_explanation)

cis<-sub_blood[sub_blood$userSet=="cis_only (107004 regions)",]
cis<-cis[cis$p.adjust<0.05&cis$oddsRatio>0,]
table(cis$seg_explanation)


trans<-sub_blood[sub_blood$userSet=="trans_only (3864 regions)",]
trans<-trans[trans$p.adjust<0.05&trans$oddsRatio>0,]
table(trans$seg_explanation)


prom<-data.frame(sub_blood[grep("Prom",sub_blood$seg_explanation),])
prom_amb<-prom[prom$userSet=="ambivalent (10941 regions)",]
prom_cis<-prom[prom$userSet=="cis_only (107004 regions)",]
prom_trans<-prom[prom$userSet=="trans_only (3864 regions)",]

prom<-data.frame(sub_blood[grep("PromD1",sub_blood$seg_explanation),])
prom[prom$p.adjust<0.05,c(1,5,22,26,27)]


####
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

