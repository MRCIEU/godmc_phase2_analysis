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

##load regiondb
regionDB0 <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")
#cistrome_cistrome  codex encode_tfbs ucsc_features cistrome_epigenome  encode_segmentation sheffield_dnase

#
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
max(clumped[which(clumped$cis==TRUE),"pval"])
#1e-4
max(clumped[which(clumped$cis==FALSE),"pval"])
#5e-8


#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(clumped$snp%in%flip[,1])
#clumped<-clumped[-w,]

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(clumped$snp%in%indels[,1]) #129
#clumped<-clumped[-w,]

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
 
#clumped<-clumped[which(clumped$cpg%in%retaincpg),]
#nrow(clumped)

clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))

#clumped<-clumped[which(clumped$pval<1e-14),]

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
  #sub[,signif:=any(p.adjust<=pval_lim),by="target"]
  #sub=sub[signif==TRUE]
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
  pl3=ggplot(sub,aes(x=target,y=-log10(p.adjust),size=logOddsRatio,fill=lineage_count,col=(cellState=="Malignant")))+geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+geom_point(alpha=0.7,shape=21,stroke=1)+facet_wrap(~userSet,scale="free_x",ncol=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+scale_size(range=c(1,4))+scale_fill_manual(values=structure(sub_colors$color, names=sub_colors$lineage_count))+scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+guides(fill = guide_legend(ncol=3,override.aes = list(size=6)))
  print(pl3)
  dev.off()
  
}

produceLOLA_plots_intern = function(grs,plot_pref,height=35,width=18,recreate=FALSE){

  #run LOLA cpgs
  simpleCache(cacheName=paste0(plot_pref,"_cpgs"),instruction="runLOLA(fg_cpg, bg_cpg, regionDB0, cores=5)",cacheDir=getwd(),recreate=recreate,assignToVariable="locResults_cpgs",buildEnvir=c(fg_cpg=grs$fg_cpg,bg_cpg=grs$bg$cpg))
   
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
data[,cpgchr:=gsub("24","Y",cpgchr),]
#data[,cpg_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]

####run with internal background
#grs_mqtl_effect=create_grs(data=data,selector=c("cpg_change=='mqtl_effect>0'","cpg_change=='mqtl_effect<0'","cpg_change=='ambivalent'","snp_change=='mqtl_effect>0'","snp_change=='mqtl_effect<0'","snp_change=='ambivalent'"))
#produceLOLA_plots_intern(grs=grs_mqtl_effect,plot_pref="mqtl_effect",height=35,width=18)


data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'"))

table(data$cpg_cis)

#ambivalent      FALSE       TRUE 
#     58803       4917     219215 

#filtered p<1e-14
#ambivalent      FALSE       TRUE 
#     51782       5421     197134 

#ambivalent      FALSE       TRUE 
#     54062       5076     210619

w<-which(data$cpg_cis=="FALSE")
dim(unique(data[w,"cpg"]))
#[1] 4260    1

w<-which(data$cpg_cis=="TRUE")
dim(unique(data[w,"cpg"]))
#[1] 157633      1

w<-which(data$cpg_cis=="ambivalent")
dim(unique(data[w,"cpg"]))
#[1] 27910     1

table(data$cis)

# FALSE   TRUE 
# 24276 258659

#FALSE   TRUE 
# 22913 246844 

#produceLOLA_plots_intern(grs=grs_cis,plot_pref="cis",height=35,width=18,recreate=TRUE)


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
m<-match(Illumina450_dt$cpgID,data$cpg,)
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

lola_res0=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)
lola_res0$logOddsRatio<-log(lola_res0$oddsRatio)
plotLOLA(locResults_all=lola_res0,plot_pref="cpg_core450kbg",height=10,width=18)

#Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
#lola_res0_matched=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB0, cores=5)
Illumina450_bg_matched=unique(with(Illumina450_bg_matched,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))

cpg_bg_gr_matched=unique(c(Illumina450_bg_matched,GoDMC_cpg_gr))

lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB0, cores=5)
lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
plotLOLA(locResults_all=lola_res0_matched,plot_pref="cpg_corebg_matched",height=10,width=18)
save(lola_res0_matched,lola_res0,file="../results/enrichments/lola_core_mqtlcpg.rdata")

##ambivalent/FALSE/TRUE

w<-which(is.na(Illumina450_dt$cpg_cis))
Illumina450_dt$cpg_cis[w]<-"0bg"
n<-names(table(Illumina450_dt$cpg_cis))[-1]

t<-table(Illumina450_dt$quantiles,Illumina450_dt$cpg_cis)
bg.matched.subset<-list()
Illumina450_bg_matchedcis<-list()
res_all<-list()
l<-list()

for (j in 1:length(n)){
cat(n[j],"\n")
w<-which(t[,(j+1)]<5)
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

bg.matched.subset[[j]]<-do.call("rbind",controllist)
Illumina450_bg_matchedcis[[j]]<-Illumina450_dt[which(Illumina450_dt$cpgID%in%bg.matched.subset[[j]]$cpgID),]
Illumina450_bg_matchedcis[[j]]<-unique(with(Illumina450_bg_matchedcis[[j]],GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
data_subset<-data[which(data$cpg_cis==n[j]),]
GoDMC_cpg_gr=unique(with(data_subset,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
l[[j]]<-length(GoDMC_cpg_gr)


cpg_bg_gr_matched=unique(c(Illumina450_bg_matchedcis[[j]],GoDMC_cpg_gr))
lola_res0_matched=runLOLA(GoDMC_cpg_gr, cpg_bg_gr_matched, regionDB0, cores=5)
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


plotLOLA(locResults_all=res_all,plot_pref="cpg_corebg_matched_cis",height=10,width=18)

save(res_all,file="../results/enrichments/lola_core_mqtlcpg_cis.rdata")


res_all2<-res_all[which(res_all$collection%in%c("codex","encode_tfbs")),]
length(unique(res_all2$antibody))
#[1] 262
res_all2$antibody<-gsub("_.*","",res_all2$antibody)
length(unique(res_all2$antibody))
#[1] 228
m<-match(res_all2$cellType,cellType_conversions$cellType)
res_all2<-data.table(data.frame(res_all2,cellType_conversions[m,]))
length(unique(res_all2$cellType_corr))
#27
length(unique(res_all2$Tissue))
#26
t1<-table(res_all2$Tissue)
dim(res_all2)
#[1] 2634   32

res_all2[which(res_all2$logOddsRatio=="-Inf"),]
table(res_all2[which(res_all2$logOddsRatio=="-Inf"),"pValueLog"])
# 0 
#10 

res_all2<-res_all2[which(res_all2$logOddsRatio!="-Inf"),]
res_all2[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
res_all2<-res_all2[which(res_all2$p.adjust<0.001),]
dim(res_all2)
#[1] 1052   33
t2<-table(res_all2$Tissue)
data.frame(t1,t2)
#             Var1 Freq         Var1.1 Freq.1
#1           Blood 1362          Blood    519
#2    Blood vessel   51   Blood vessel     29
#3            Bone    9           Bone      5
#4           Brain   78          Brain     37
#5      Cerebellum    3     Cerebellum      2
#6          Cervix  195         Cervix     59
#7  Choroid plexus    3 Choroid plexus      2
#8           Colon   21          Colon     11
#9        Endoderm    6       Endoderm      3
#10    Endometrium   24    Endometrium      6
#11            ESC  240            ESC     86
#12      Esophagus    3      Esophagus      2
#13            Eye    3            Eye      2
#14        Gingiva    3        Gingiva      2
#15          Heart    6          Heart      4
#16         Kidney   36         Kidney     15
#17          Liver  228          Liver     89
#18           Lung  141           Lung     68
#19          Mamma  129          Mamma     68
#20         Muscle   12         Muscle      4
#21       Pancreas   12       Pancreas      5
#22       Placenta    3       Placenta      2
#23         Retina    3         Retina      2
#24           Skin   51           Skin     27
#25    Spinal cord    3    Spinal cord      2
#26         Testis    9         Testis      1

load("tfbsbycpgmean.Robj")
load("../results/enrichments/lola_core_mqtlcpg.rdata")
res_all<-lola_res0_matched
res_all$antibody<-gsub("_.*","",res_all$antibody)

df.out$filename<-gsub(".bed","",df.out$antibody)
w<-which(names(df.out)%in%"antibody")
df.out<-df.out[,-w]
#df.out$antibody2<-gsub("\\_.*","",df.out$antibody)
#df.out$cellType<-sub(".*._", "", df.out$antibody)

w<-which(df.out$filename%in%res_all$filename)
df.out2<-df.out[w,]
m<-match(df.out2$filename,res_all$filename)
df.out2<-data.frame(df.out2,res_all[m,])

df.out0<-df.out2[which(df.out2[,1]==0),]
df.out1<-df.out2[which(df.out2[,1]==1),]
df.out1<-data.table(df.out1)
##  
pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=df.out1,cellType_conversions=cellType_conversions)
  
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
  

pdf(paste0("meanmethylation_tfbs_cis.pdf"),height=10,width=18)
  p1<-ggplot(sub,aes(x=antibody,y=mean,fill=lineage_count,col=(cellState=="Malignant")))+
  geom_hline(yintercept=mean(df.out0$mean),col="black",linetype="dashed")+
  geom_point(size=4,alpha=0.7,shape=21,stroke=1)+
  #facet_wrap(~userSet,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+scale_fill_manual(values=structure(sub_colors$color, names=sub_colors$lineage_count))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))+guides(fill = guide_legend(ncol=3))
  labs(y="mean methylation")
  print(p1)
  dev.off()


m<-match(df.out1$cellType,cellType_conversions$cellType)
df.out1<-data.table(data.frame(df.out1,cellType_conversions[m,]))

hypo<-df.out1[df.out1$mean<0.1,]
length(unique(hypo$cellType_corr))
#18
length(unique(hypo$Tissue))
#15
t1<-table(hypo$Tissue)
dim(hypo)
#[1] 2634   32d
df<-unique(data.frame(hypo$antibody,hypo$Tissue))
df<-data.frame(table(df[,1]))
dim(df) #127
dim(df[df$Freq>5,]) 
hypo[which(hypo$antibody%in%c("Nrf1")),c("mean","Lineage","Tissue","cellType")]
hypo[which(hypo$antibody%in%c("BRCA1")),c("mean","Lineage","Tissue","cellType")]

#nvars<-names(ov.all)[6:ncol(ov.all)]
#ov.all2<-ov.all[,seq(2,6:ncol(ov.all)),with=FALSE]
#row.names(ov.all2)<-ov.all[,1]
#out<-ov.all%>%group_by(paste0(nvars))%>%summarise(mean=mean(meancpg))
#out<-ov.all[, lapply(.SD,mean(meancpg)), by=ov.all[,6]]
#out<-ov.all[, lapply(.SD,mean), by=eval(colnames(ov.all)[6:ncol(ov.all)])]
#.SD is the (S)ubset of (D)ata excluding group columns



res_all2[userSet=="trans_only (3864 regions)",][collection=="encode_tfbs"|collection=="codex", ][1:15,]
res_all2[userSet=="ambivalent (10941 regions)",][collection=="encode_tfbs"|collection=="codex", ][1:15,]
res_all2[userSet=="ambivalent (10941 regions)",][antibody=="CTCF"][collection=="encode_tfbs"|collection=="codex", ][1:15,]
res_all2[userSet=="ambivalent (10941 regions)",][antibody=="CTCF",][collection=="encode_tfbs"|collection=="codex","cellType_corr"] 
dim(unique(res_all2[userSet=="ambivalent (10941 regions)",][antibody=="CTCF",][collection=="encode_tfbs"|collection=="codex","cellType_corr"] ))
#19
dim(unique(res_all2[userSet=="ambivalent (10941 regions)",][antibody=="Rad21",][collection=="encode_tfbs"|collection=="codex","cellType_corr"] ))
#7
dim(unique(res_all2[userSet=="ambivalent (10941 regions)",][antibody=="SMC3",][collection=="encode_tfbs"|collection=="codex","cellType_corr"] ))
#5

amb<-unique(res_all2[userSet=="ambivalent (10941 regions)",][collection=="encode_tfbs"|collection=="codex"])
df<-unique(data.frame(amb$antibody,amb$cellType_corr))
df<-data.frame(table(df[,1]))
dim(df) #127
dim(df[df$Freq>3,]) #12

#df[df$Freq>5,]
#    Var1 Freq
#15 CEBPB    6
#21  CTCF   19
#92 Rad21    7



df<-unique(data.frame(amb$antibody,amb$Tissue))
df<-data.frame(table(df[,1]))
dim(df) #126
dim(df[df$Freq>3,]) #11

tf<-df[df$Freq>3,"Var1"]
for (i in 1:length(tf)){
lineage<-unique(res_all2[userSet=="ambivalent (10941 regions)",][antibody==tf[i],][collection=="encode_tfbs"|collection=="codex",c("Lineage")])
lineage<-data.frame(lineage)
w<-which(lineage[,1]%in%c("Endodermal","Ectodermal","Mesodermal"))
if(length(w)==3){
print(tf[i])  
}
}

df<-data.frame(table(res_all2$Tissue))
#516/1038
or<-ddply(amb,~antibody,summarise,mean=mean(oddsRatio))
mean(or$mean)
#1.335355

trans<-unique(res_all2[userSet=="trans_only (3864 regions)",][collection=="encode_tfbs"|collection=="codex"])
df<-unique(data.frame(trans$antibody,trans$cellType_corr))
df<-data.frame(table(df[,1]))
dim(df) #168
dim(df[df$Freq>3,]) #36
or<-ddply(trans,~antibody,summarise,mean=mean(oddsRatio))
mean(or$mean)
#[1]  1.54232

df<-unique(data.frame(trans$antibody,trans$Tissue))
df<-data.frame(table(df[,1]))
dim(df) #168
dim(df[df$Freq>3,]) #19

#df[which(df$Freq>5),]
#         Var1 Freq
#31      c-Myc    6
#34       CTCF   23
#115      Pol2   12
#116  Pol2-4H8    7
#123     Rad21    6
#133 Sin3Ak-20    6
#145      TAF1    6
#164       YY1    7

tf<-df[df$Freq>3,"Var1"]
for (i in 1:length(tf)){
lineage<-unique(res_all2[userSet=="ambivalent (10941 regions)",][antibody==tf[i],][collection=="encode_tfbs"|collection=="codex",c("Lineage")])
lineage<-data.frame(lineage)
w<-which(lineage[,1]%in%c("Endodermal","Ectodermal","Mesodermal"))
if(length(w)==3){
print(tf[i])  
}
}

cis<-unique(res_all2[userSet=="cis_only (107004 regions)",][collection=="encode_tfbs"|collection=="codex"])
df<-unique(data.frame(cis$antibody,cis$cellType_corr))
df<-data.frame(table(df[,1]))
dim(df) #17
dim(df[df$Freq>3,]) #0

df<-unique(data.frame(cis$antibody,cis$Tissue))
df<-data.frame(table(df[,1]))
dim(df) #17
dim(df[df$Freq>3,]) #0

or<-ddply(cis,~antibody,summarise,mean=mean(oddsRatio))
mean(or$mean)
#1.35719

###

###
#


