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
library(ggrepel)

theme_set(theme_bw())

regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

##load cell type conversion and colors
cellType_conversions=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/CellTypes.tsv",drop="collection")
colors=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/color.tsv")

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
  #LOLA_res[description=="T-cell acute lymphoblastic leukaemia (T-ALL) cell line.",c("Lineage1","Lineage","cellType_corr"):=list(Lineage1="Lymphoid",Lineage="Lymphoid",cellType_corr="T lymphocyte"),]
  #LOLA_res[,lineage_count_allstate:=length(unique(filename[!is.na(filename)])),by=c("Lineage","cellState")]
  #LOLA_res[,lineage_count_all:=length(unique(filename[!is.na(filename)])),by=c("Lineage")]
  
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
  
  pdf(paste0(plot_pref,"_OR_All.pdf"),height=10,width=16)
  m<-max(sub$oddsRatio)+1
  pl4<-ggplot(res_all2,aes(x=antibody,y=oddsRatio,size=pValueLog))+
  geom_hline(yintercept=1,col="black",linetype="dotted")+
  geom_point(aes(color=Tissue))+
    facet_wrap(~userSet,scale="free_x",ncol=1)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
    scale_size(range=c(1,4))+
    guides(fill = guide_legend(ncol=20))+
    scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
    labs(y="Odds ratio (log scale)",x="Transcription factor") +
    scale_fill_brewer(type="qual") +
    theme(legend.text=element_text(size=12))
  
    print(pl4)
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

####run with external background for CpGs
m<-match(Illumina450_dt$cpgID,data$cpg)
df<-data[m,"cluster"]
Illumina450_dt$cluster<-df

#Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
#Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)

Illumina450_dt$GC_freqquantile<-cut(Illumina450_dt$GC_freq, breaks=c(quantile(Illumina450_dt$GC_freq,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)
Illumina450_dt$CpG_freqquantile<-cut(Illumina450_dt$CpG_freq, breaks=c(quantile(Illumina450_dt$CpG_freq,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)

Illumina450_dt$quantiles<-paste(Illumina450_dt$CpG_freqquantile,Illumina450_dt$GC_freqquantile)

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
if(t[i,1]>5&t[i,(j+1)]>0){
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
lola_res0_matched$userSet<-n[j]
res_all<-rbind(res_all,lola_res0_matched)

}
save(res_all,file="tfbs_enrichment.matched.Robj")
###DONE UNTIL HERE

#plotLOLA(locResults_all=res_all,plot_pref="cpg_corebg_matched_cis_updated",height=10,width=18)


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

w1<-which(res_all2$logOddsRatio=="-Inf")
table(res_all2$pValueLog[w1])

w2<-which(res_all2$logOddsRatio=="Inf")
table(res_all2$pValueLog[w2])

res_all2<-res_all2[which(res_all2$logOddsRatio!="-Inf"),]
res_all2[,Pvalue:=10^(-pValueLog)]
res_all2[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
res_all2<-data.frame(res_all2)

locResults_all=res_all

plot_pref="cpg_corebg_matched_cis_updated"
pval_lim=0.001
  
  locResults=process_LOLA(LOLA_res=data.table(res_all2),cellType_conversions=cellType_conversions)
   

  sub=locResults
  sub[,signif:=any(p.adjust<=pval_lim),by="target"]
  sub=sub[signif==TRUE]
  sub[which(sub$logOddsRatio!="Inf"),]
sub<-sub[which(sub$p.adjust<0.001),]

res_all2<-res_all2[which(res_all2$p.adjust<0.001),]

res_all3<-res_all2[,c("userSet","antibody","oddsRatio","Pvalue","p.adjust","support","b","c","d","cellType","cellType_corr","Tissue")]
names(res_all3)<-c("Community","Type","OR","Pvalue","FDR_Pvalue","support","b","c","d","CellType","cellType_category","Tissue")                                                   
write.table(res_all3,"TableSXX_LOLA_tfbs_communities.txt",sep="\t",quote=F,row.names=F,col.names=T)

#pdf(paste0(plot_pref,"_OR_All.pdf"),height=10,width=18)
locResults_sig <- subset(locResults, p.adjust < 0.001)
p<-unique(paste(locResults_sig$userSet,locResults_sig$target))
m <-match(p,paste(locResults_sig$userSet,locResults_sig$target))
lab<-locResults_sig[m,c("target","userSet","p.adjust","oddsRatio")]

  w1<-which(locResults$oddsRatio=="Inf")
  locResults2<-locResults[-w1,]
  m1<-max(locResults2$oddsRatio)+1
  m2<-min(locResults2$oddsRatio)-1
  pl4<-ggplot(locResults2,aes(x=target,y=oddsRatio,size=pValueLog,label=userSet))+
    geom_hline(yintercept=1,col="black",linetype="dotted")+
    geom_point(aes(colour=factor(userSet),label=userSet))+
    guides(colour = guide_legend(nrow = 3)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
    #scale_size(range=c(1,4))+
    geom_label_repel(data=lab, aes(x=target,y=oddsRatio,label=userSet), size=4) +
    scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20,50,100)) +
    theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Transcription factor", fill="community")
    pl4$labels$colour<-"Community"
   
    ggsave(pl4, file=paste0(plot_pref,"_OR_All.pdf"),height=10,width=18)
 # print(pl4)
 # dev.off()

res16<-subset(res_all,userSet==16)
res16[,Pvalue:=10^(-pValueLog)]
res16[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
res16[res16$Pvalue<0.05,]
