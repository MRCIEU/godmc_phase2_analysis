library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(biovizBase)

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())


setwd("/scratch/lab_bock/jklughammer/projects/otherProjects/GoDMC/")

##load cell type conversion and colors
cellType_conversions=fread("/data/groups/lab_bock/jklughammer/gitRepos/otherProjects/GoDMC/LOLA_annot/CellTypes.tsv",drop="collection")
colors=fread("/data/groups/lab_bock/jklughammer/gitRepos/otherProjects/GoDMC/LOLA_annot/color.tsv")

##load regiondb
regionDB = loadRegionDB("/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/")

##load GoDMC data
load("for_christoph.rdata")
data=as.data.table(dat)
load("snpsetsforJohanna170725.RData")
snp_bg=as.data.table(controlset)


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
  fg_snp_list=GRangesList()
  for (sel in selector){
    fg_snp=with(unique(data[eval(parse(text=sel)),c("snpchr","ld80_start","ld80_end"),with=FALSE]),GRanges(seqnames=snpchr,ranges=IRanges(ld80_start,ld80_end))) 
    fg_snp_list[[paste0(sel,"_snp ",length(fg_snp)," regons")]]= fg_snp
  }
  
  bg_list=GRangesList()
  bg_list[["cpg"]]=with(unique(data[,c("cpgchr","cpgstart","cpgend"),with=FALSE]),GRanges(seqnames=cpgchr,ranges=IRanges(cpgstart,cpgend)))
  bg_list[["snp"]]=with(unique(data[,c("snpchr","ld80_start","ld80_end"),with=FALSE]),GRanges(seqnames=snpchr,ranges=IRanges(ld80_start,ld80_end)))
  
  return(list(fg_cpg=fg_cpg_list,fg_snp=fg_snp_list,bg=bg_list))
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
  #run LOLA snps
  simpleCache(cacheName=paste0(plot_pref,"_snps"),instruction="runLOLA(fg_snp, bg_snp, regionDB, cores=5)",cacheDir=getwd(),recreate=recreate,assignToVariable="locResults_snps",buildEnvir=c(fg_snp=grs$fg_snp,bg_snp=grs$bg$snp))
  
  #combine and process
  locResults_all=rbindlist(list(locResults_cpgs,locResults_snps))
  plotLOLA(locResults_all=locResults_all,plot_pref=plot_pref,height=height,width=width)
}

##prepare GoDMC data: merge CpGs and SNPs that are in proximity to eachother to avoid infalting the results, 1kb around cp
Illumina450_sub=Illumina450_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina450_sub,c("cpgID","pos"),c("cpgname","ill_pos"))
data=merge(data,Illumina450_sub,by="cpgname",all.x=TRUE)
#check pos (should be all true)
table(data[,ill_pos==cpgpos,])

data[,cpgchr:=gsub("23","X",cpgchr),]
data[,snpchr:=gsub("23","X",cpgchr),]


data[,cpg_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
data[,snp_change:=ifelse(all(mqtl_effect>0),"mqtl_effect>0",ifelse(all(mqtl_effect<0),"mqtl_effect<0","ambivalent")),by=c("snpchr","ld80_start","ld80_end")]

data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgstart","cpgend")]
data[,snp_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("snpchr","ld80_start","ld80_end")]


####run with internal background
grs_mqtl_effect=create_grs(data=data,selector=c("cpg_change=='mqtl_effect>0'","cpg_change=='mqtl_effect<0'","cpg_change=='ambivalent'","snp_change=='mqtl_effect>0'","snp_change=='mqtl_effect<0'","snp_change=='ambivalent'"))
produceLOLA_plots_intern(grs=grs_mqtl_effect,plot_pref="mqtl_effect",height=35,width=18)

grs_cis=create_grs(data=data,selector=c("cpg_cis=='TRUE'","cpg_cis=='FALSE'","cpg_cis=='ambivalent'","snp_cis=='TRUE'","snp_cis=='FALSE'","snp_cis=='ambivalent'"))
produceLOLA_plots_intern(grs=grs_cis,plot_pref="cis",height=35,width=18,recreate=TRUE)

####run with external background for CpGs

hg19_Illumina450_gr=with(Illumina450_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
seq_Illumina450=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina450_gr)
#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~3000)
Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

Illumina450_dt[,GC_freq:=letterFrequency(seq_Illumina450, "CG", as.prob=T),]
Illumina450_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina450, step=2, as.prob=T)[,"CG"],]
Illumina450_dt[,isGoDMC:=ifelse(cpgID%in%data$cpgname,TRUE,FALSE),]

#plot CG and CpG frequency for GoDMC cpgs and background 
pdf("compare_seq_properties.pdf",height=3,width=4)
ggplot(Illumina450_dt,aes(x=GC_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq,col=isGoDMC))+geom_density()
ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()
dev.off()


GoDMC_cpg_gr=unique(with(data,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
#use all for background
Illumina450_bg=unique(with(Illumina450_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
lola_res=runLOLA(GoDMC_cpg_gr, Illumina450_bg, regionDB, cores=5)
plotLOLA(locResults_all=lola_res,plot_pref="cpg_extbg",height=10,width=18)

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

####run with external background for SNPs


GoDMC_snp_gr=unique(with(data,GRanges(seqnames = Rle(snpchr), IRanges(start=ld80_start, end=ld80_end),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()
for (i in 1:10){
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

for (i in 1:10){
lola_res=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], regionDB, cores=5)
plotLOLA(locResults_all=lola_res,plot_pref=paste0("snp_extbg_",i),height=6,width=9)
}



####annotate CpGs and SNPs with LOLA regions
getEnrichmentOverlaps=function(EOI,userset,regionDB){
  enriched_regions=GRangesList()
  avail_names=vector()
  for (i in 1:nrow(EOI)){
    extr=extractEnrichmentOverlaps(EOI[i,], userset, regionDB)
    extr$target=EOI[i,]$target
    enriched_regions=append(enriched_regions,GRangesList(extr))
    avail_names=c(avail_names,EOI[i,]$filename)
  }
  names(enriched_regions)=avail_names
  return(enriched_regions)
}

annotate_gr = function(gr,collections=c("encode_tfbs","codex","cistrome_epigenome")){
  all_res=runLOLA(gr, gr, regionDB, cores=5)
  EOI=all_res[support>0&collection%in%collections]
  enriched_regions=getEnrichmentOverlaps(EOI,gr,regionDB)
  
  enriched_regions_dt=as.data.table(enriched_regions)
  setnames(enriched_regions_dt,"group_name","filename")
  enriched_regions_dt=merge(enriched_regions_dt,EOI[,c("filename","cellType","antibody"),with=FALSE],by="filename")
  
  return(enriched_regions_dt)
}


all_cpg_gr=unique(with(data,GRanges(seqnames = Rle(cpgchr), IRanges(start=cpgpos, end=cpgpos),strand=Rle("*"),ID=cpgname)))
all_snp_gr=unique(with(data,GRanges(seqnames = Rle(cpgchr), IRanges(start=snppos, end=snppos),strand=Rle("*"),ID=snpname)))


cpgs_annot=annotate_gr(all_cpg_gr)
snps_annot=annotate_gr(all_snp_gr)

write.table(cpgs_annot,"annotated_cpgs.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(snps_annot,"annotated_snps.tsv",sep="\t",quote=FALSE,row.names=FALSE)

save(cpgs_annot,file="annotated_cpgs.RData")
save(snps_annot,file="annotated_snps.RData")


#back-check some of the annotations by selecting different lines (result should be one line) --> seems ok
check_sub=snps_annot[8888,]
filename=as.character(check_sub$filename)
cseqnames=as.character(check_sub$seqnames)
cstart=as.numeric(check_sub$start)
cend=as.numeric(check_sub$end)
as.data.table(getRegionSet(regionDB=regionDB,filenames=filename))[seqnames==cseqnames&start<cstart&end>cend]
######################




##CpG foreground: 1kb around CpG
##CpG background: 450k pos, 1kb, matched CG, CpG, and rep

##SNP foreground: LD80 regions
##SNP background: positive vs negative effectsize

#analysis: cis-trans, Isq (meta analysis heterogeneity),pchic (celltypes of chromosome interaction between mQTL and SNP)
