library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(biovizBase)
library(meffil)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

##load cell type conversion and colors
cellType_conversions=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/CellTypes.tsv",drop="collection")
colors=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/color.tsv")

setwd("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments")

##load cell type conversion and colors
cellType_conversions=fread("/data/groups/lab_bock/jklughammer/gitRepos/otherProjects/GoDMC/LOLA_annot/CellTypes.tsv",drop="collection")
colors=fread("/data/groups/lab_bock/jklughammer/gitRepos/otherProjects/GoDMC/LOLA_annot/color.tsv")

##load regiondb
regionDB0 <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
max(clumped[which(clumped$cis==TRUE),"pval"])
#1e-4
max(clumped[which(clumped$cis==FALSE),"pval"])
#5e-8

flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(clumped$snp%in%flip[,1])
clumped<-clumped[-w,]

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
 
clumped<-clumped[which(clumped$cpg%in%retaincpg),]
nrow(clumped)

#head(dat)
#  cpgchr    cpgpos    cpgname snpchr    snppos              snpname
#1   chr1 230560793 cg00000363   chr1 230638673   chr1:230638673:SNP
#2   chr1 166958439 cg00001349   chr1 166953284   chr1:166953284:SNP
#3   chr1 166958439 cg00001349   chr1 166319590 chr1:166319590:INDEL
#4   chr1 200011786 cg00001583   chr1 199963470   chr1:199963470:SNP
#5   chr1 153515502 cg00002646   chr1 153484034   chr1:153484034:SNP
#6   chr1 169396706 cg00002719   chr1 169358424   chr1:169358424:SNP
#  effect_allele other_allele effect_allele_freq mqtl_effect mqtl_pval
#1             a            c             0.0612   0.2209044 1.025e-07
#2             t            c             0.2456  -0.2167937 2.828e-27
#3             d            i             0.0195   0.6823649 4.763e-06
#4             a            g             0.8575   0.2051145 1.632e-16
#5             t            c             0.8907  -0.1466724 6.052e-07
#6             a            g             0.1224  -0.2381370 1.886e-19
#  meta_directions  Isq samplesize  cis ld80_start  ld80_end ld80_proxies
#1         ++++??+  0.0       5400 TRUE  230638673 230638673            0
#2         ------- 66.4       6673 TRUE  166953455 166958937           14
#3         ?+?????  0.0       1177 TRUE  166327630 166327630            1
#4         +++++++  0.0       6673 TRUE  199963622 199967877            4
#5         ---?---  0.0       5954 TRUE  153484148 153512347           44
#6         -------  0.0       6673 TRUE  169361319 169447279           30
#  ld80_dist                            code
#1         0   cg00000363 chr1:230638673:SNP
#2      5482   cg00001349 chr1:166953284:SNP
#3         0 cg00001349 chr1:166319590:INDEL
#4      4255   cg00001583 chr1:199963470:SNP
#5     28199   cg00002646 chr1:153484034:SNP
#6     85960   cg00002719 chr1:169358424:SNP
#                                                                           pchic
#1                                                                           <NA>
#2                                                                           <NA>
#3                                                                           <NA>
#4                                                                           <NA>
#5 Mon,Mac0,Mac1,Mac2,Neu,MK,MK,EP,EP,Ery,Ery,FoeT,nCD4,tCD4,aCD4,naCD4,nCD8,tCD8
#6                                                                           <NA>




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

process_LOLA = function (LOLA_res, collections=c("codex","encode_tfbs"),cellType_conversions){
  
  LOLA_res=LOLA_res[!is.na(userSet)]
  
  LOLA_res=LOLA_res[collection%in%collections]
  #changed form exp() to 10^ to accomodate change in LOLA  	
  LOLA_res[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
  LOLA_res[,mlog10p.adjust:=-log10(p.adjust),]
  
  ##standardize cellTypes
  LOLA_res[cellType==""|is.na(cellType),cellType:="Not defined",]
  LOLA_res=merge(LOLA_res,cellType_conversions,by="cellType",all.x=TRUE)
  ##correct wrong annotation
  LOLA_res[description=="T-cell acute lymphoblastic leukaemia (T-ALL) cell line.",c("Lineage1","Lineage","cellType_corr"):=list(Lineage1="Lymphoid",Lineage="Lymphoid",cellType_corr="T lymphocyte"),]
  LOLA_res[,lineage_count_allstate:=length(unique(filename[!is.na(filename)])),by=c("Lineage","cellState")]
  LOLA_res[,lineage_count_all:=length(unique(filename[!is.na(filename)])),by=c("Lineage")]
  
  #standardize antibodies
  LOLA_res[,target:=toupper(sub("-","",unlist(lapply(antibody,function(x){spl=unlist(strsplit(x,"_|eGFP-"));spl[spl!=""][1]})))),]
  
  return(LOLA_res)  
  
}

plotLOLA=function(locResults_all,plot_pref,height=35,width=18){
  
  ##process LOLA results
  pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=locResults_all,cellType_conversions=cellType_conversions)
  
  sub=locResults
  sub[,signif:=any(p.adjust<=pval_lim),by="target"]
  sub=sub[signif==TRUE]
  sub[,lineage_count:=paste0(Lineage," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",lineage_count_all),by=c("Lineage")]
  sub[,lineage_count_state:=paste0(Lineage," ",length(unique(filename[which(p.adjust<=pval_lim)])),"/",lineage_count_allstate),by=c("Lineage","cellState")]
  
  sub_colors=sub[,.N,by=c("Lineage","lineage_count")]
  sub_colors=merge(sub_colors,colors,by="Lineage")
  
  sub_colors_state=sub[,.N,by=c("Lineage","cellState","lineage_count_state")]
  sub_colors_state=merge(sub_colors_state,colors,by="Lineage")
  
  ##summary plot
  sub_sum=sub[,list(mean_mlog10p.adjust=mean(mlog10p.adjust[p.adjust<=pval_lim&mlog10p.adjust!=Inf]),count=length(userSet[p.adjust<=pval_lim])),by="userSet"]
  
  pdf(paste0(plot_pref,"_summary.pdf"),height=height/7,width=10)
  pl=ggplot(sub_sum,aes(x=userSet,y=count))+geom_bar(stat="identity",aes(fill=mean_mlog10p.adjust))+ylab("Enriched region sets")+coord_flip()+scale_fill_gradient(low="blue",high="red")
  print(pl)
  dev.off()
  
  ##bubble plots
  pdf(paste0(plot_pref,"_Malignant.pdf"),height=height,width=width)
  pl1=ggplot(sub[cellState=="Malignant"],aes(x=target,y=-log10(p.adjust),size=logOddsRatio,fill=lineage_count_state))+geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+geom_point(alpha=0.7,shape=21,stroke=1)+facet_wrap(~userSet,scale="free_x",ncol=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+scale_size(range=c(1,4))+scale_fill_manual(values=structure(sub_colors_state[cellState=="Malignant"]$color, names=sub_colors_state[cellState=="Malignant"]$lineage_count))+guides(fill = guide_legend(ncol=3))
  print(pl1)
  dev.off()
  
  pdf(paste0(plot_pref,"_Healthy.pdf"),height=height,width=width/3)
  pl2=ggplot(sub[cellState=="Healthy"],aes(x=target,y=-log10(p.adjust),size=logOddsRatio,fill=lineage_count_state))+geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+geom_point(alpha=0.7,shape=21,stroke=1)+facet_wrap(~userSet,scale="free_x",ncol=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+scale_size(range=c(1,4))+scale_fill_manual(values=structure(sub_colors_state[cellState=="Healthy"]$color, names=sub_colors_state[cellState=="Healthy"]$lineage_count))+guides(fill = guide_legend(ncol=3))
  print(pl2)
  dev.off()
  
  pdf(paste0(plot_pref,"_All.pdf"),height=height,width=width+2)
  pl3=ggplot(sub,aes(x=target,y=-log10(p.adjust),size=logOddsRatio,fill=lineage_count,col=(cellState=="Malignant")))+geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+geom_point(alpha=0.7,shape=21,stroke=1)+facet_wrap(~userSet,scale="free_x",ncol=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+scale_size(range=c(1,4))+scale_fill_manual(values=structure(sub_colors$color, names=sub_colors$lineage_count))+scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+guides(fill = guide_legend(ncol=3))
  print(pl3)
  dev.off()
  
}

produceLOLA_plots_intern = function(grs,plot_pref,height=35,width=18,recreate=FALSE){

  #run LOLA cpgs
  simpleCache(cacheName=paste0(plot_pref,"_cpgs"),instruction="runLOLA(fg_cpg, bg_cpg, regionDB, cores=5)",cacheDir=getwd(),recreate=recreate,assignToVariable="locResults_cpgs",buildEnvir=c(fg_cpg=grs$fg_cpg,bg_cpg=grs$bg$cpg))
   
  #combine and process
  locResults_all=locResults_cpgs
  locResults_all$logOddsRatio<-log(locResults_all$oddsRatio)
  plotLOLA(locResults_all=locResults_all,plot_pref=plot_pref,height=height,width=width)
}

##prepare GoDMC data: merge CpGs and SNPs that are in proximity to eachother to avoid infalting the results, 1kb around cp
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpg","ill_pos"))
data=merge(data,Illumina450_sub,by="cpg",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==cpgpos,])

data[,cpgchr:=gsub("23","X",cpgchr),]
#data[,cpg_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]

####run with internal background
grs_mqtl_effect=create_grs(data=data,selector=c("cpg_change=='mqtl_effect>0'","cpg_change=='mqtl_effect<0'","cpg_change=='ambivalent'","snp_change=='mqtl_effect>0'","snp_change=='mqtl_effect<0'","snp_change=='ambivalent'"))
produceLOLA_plots_intern(grs=grs_mqtl_effect,plot_pref="mqtl_effect",height=35,width=18)

data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'"))
produceLOLA_plots_intern(grs=grs_cis,plot_pref="cis",height=35,width=18,recreate=TRUE)


####run with external background for CpGs

hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)
#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~3000)
Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isGoDMC:=ifelse(cpgID%in%data$cpg,TRUE,FALSE),]
Illumina450_dt<-Illumina450_dt[Illumina450_dt$cpgID%in%retaincpg,]

Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
Illumina450_dt$quantiles<-paste(Illumina450_dt$CpG_freqquantile,Illumina450_dt$GC_freqquantile)

t<-table(Illumina450_dt$quantiles,Illumina450_dt$isGoDMC)
w<-which(t[,2]<5)
pr<-data.frame(t[-w,1]/t[-w,2])
m<-min(pr)*t[,2]

controllist<-list()
for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[i,1]>5&t[i,2]>0){
subgroup<-Illumina450_dt[which(Illumina450_dt$quantiles==row.names(t)[i]&Illumina450_dt$isGoDMC==FALSE),]
id<-subgroup[sample(nrow(subgroup), size=round(m[i],0), replace=FALSE),"cpgID"]
controllist[[i]]<-id
}}

bg.matched<-do.call("rbind",controllist)
Illumina450_bg_matched<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched$cpgID),]
Illumina450_godmc<-Illumina450_dt[which(Illumina450_dt$cpgID%in%data$cpg),]
bg_matched<-rbind(Illumina450_bg_matched,Illumina450_godmc)

GoDMC_cpg_gr=unique(with(data,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

#plot CG and CpG frequency for GoDMC cpgs and background 
pdf("compare_seq_properties.pdf",height=3,width=4)
ggplot(Illumina450_dt,aes(x=GC_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()
dev.off()

pdf("compare_seq_properties_matched.pdf",height=3,width=4)
ggplot(bg_matched,aes(x=GC_freq,col=isGoDMC))+geom_density()
ggplot(bg_matched,aes(x=GC_freq))+geom_density()
ggplot(bg_matched,aes(x=CpG_freq,col=isGoDMC))+geom_density()
ggplot(bg_matched,aes(x=CpG_freq))+geom_density()
dev.off()


#use all for background
Illumina450_bg=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
lola_res0=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)
lola_res0$logOddsRatio<-log(lola_res0$oddsRatio)

plotLOLA(locResults_all=lola_res0,plot_pref="cpg_core450kbg",height=10,width=18)

#Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
#lola_res0_matched=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)

cpg_bg_gr_matched=unique(c(Illumina450_bg_matched,GoDMC_cpg_gr))
lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB0, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
plotLOLA(locResults_all=lola_res0_matched,plot_pref="cpg_corebg_matched",height=10,width=18)

###
regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")

lola_res=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB, cores=5)

lola_res<-rbind(lola_res0,lola_res)
lola_res$logOddsRatio<-log(lola_res$oddsRatio)
plotLOLA(locResults_all=lola_res,plot_pref="cpg_extbg",height=10,width=18)
lola_res[userSet==1,][order(maxRnk, decreasing=FALSE),][1:10,]

###
#


#use gc and CpG matched background --> not really possible because half of all cpgs are GoDMC CpGs
Illumina450_dt_forbg=Illumina450_dt[order(GC_freq, CpG_freq)]
Illumina450_dt_forbg[,rank:=1:nrow(Illumina450_dt_forbg),]
bg_rank=Illumina450_dt_forbg[isGoDMC==FALSE]$rank
fg_rank=Illumina450_dt_forbg[isGoDMC==TRUE]$rank

#This takes really long
for (i in c(1:nrow(Illumina450_dt_forbg))){
  if (Illumina450_dt_forbg[i]$isGoDMC==TRUE){
    min_idx=which.min(abs(bg_rank-Illumina450_dt_forbg[i]$rank))
    Illumina450_dt_forbg[i,matched_index:=bg_rank[min_idx]]
    bg_rank=bg_rank[min_idx:length(bg_rank)]
  }else{
    Illumina450_dt_forbg[i,matched_index:=as.integer(0)]
  } 
}
