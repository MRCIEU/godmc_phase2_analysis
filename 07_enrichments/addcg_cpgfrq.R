arguments<-commandArgs(T)
no<-as.numeric(arguments[1])

library(tidyverse)
library(dplyr)
library(meffil)
library(data.table)
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicFeatures")

load("../results/enrichments/snpcontrolsets.rdata")

r.all<-data.frame()

for (i in (1:23)){
cat(i,"\n")
p<-paste0("chr",i)
r<-f.all[f.all$snpchr==p,]
r$strand<-"+"
r$snpchr<-gsub("chr23","chrX",r$snpchr)
r_dt=as.data.table(r)


#r_dt[,snpstart_pre:=ifelse(strand=="-",snppos-500,snppos-499),]
#r_dt[,snpend_pre:=ifelse(strand=="-",snppos+500,snppos+501),]

#collapse overlaps

gr_range = with(r_dt,GRanges(seqnames=snpchr,ranges=IRanges(min,max)))
gr_snp = with(r_dt,GRanges(seqnames=snpchr,ranges=IRanges(snppos,snppos)))
overlap=as.data.table(findOverlaps(gr_snp, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

r_dt[,snpstart:=start(gr_range[overlap_red$subjectHit])]
r_dt[,snpend:=end(gr_range[overlap_red$subjectHit])]
r_dt[,NsubjectHits:=overlap_red$NsubjectHits]

hg19_gr=with(r_dt, GRanges(seqnames = Rle(snpchr), IRanges(start=snpstart, end=snpend),strand=Rle(strand),ID=snppos))

seq_hg19=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_gr)

#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~ 3000 )
#Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

r_dt[,GC_freq:=letterFrequency(seq_hg19, "CG", as.prob=T),]
r_dt[,CpG_freq:=dinucleotideFrequency(seq_hg19, step=2, as.prob=T)[,"CG"],]

r<-data.frame(r_dt)
r.all<-rbind(r.all,r)
}
save(r.all,file=paste("../results/enrichments/snpcontrolsetsGC_CpGcontent.rdata",sep=""))


