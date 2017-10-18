library(data.table)
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicFeatures")


bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
names(bim)<-c("chr","snp","cm","pos","a1","a2")
p<-paste("chr",bim$chr,sep="")
p<-gsub("chr23","chrX",p)
bim$chr<-p

bim$strand<-"+"
bim_dt=as.data.table(bim)
bim_dt[,snpstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
bim_dt[,snpend_pre:=ifelse(strand=="-",pos+500,pos+501),]

#get 450k locations
#Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
#Illumina450_dt=as.data.table(Illumina450)
#Illumina450_dt[,cpgID:=row.names(Illumina450),]
#Illumina450_dt[,cpgstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
#Illumina450_dt[,cpgend_pre:=ifelse(strand=="-",pos+500,pos+501),]

#collapse overlaps
gr_range = with(bim_dt,GRanges(seqnames=chr,ranges=IRanges(snpstart_pre,snpend_pre)))
gr_cpg = with(bim_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

bim_dt[,snpstart:=start(gr_range[overlap_red$subjectHit])]
bim_dt[,snpend:=end(gr_range[overlap_red$subjectHit])]
bim_dt[,NsubjectHits:=overlap_red$NsubjectHits]

hg19_gr=with(bim_dt, GRanges(seqnames = Rle(chr), IRanges(start=snpstart, end=snpend),strand=Rle(strand),ID=snp))

seq_hg19=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_gr)

#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~ 3000 )
#Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

bim_dt[,GC_freq:=letterFrequency(seq_hg19, "CG", as.prob=T),]
bim_dt[,CpG_freq:=dinucleotideFrequency(seq_hg19, step=2, as.prob=T)[,"CG"],]

load("../results/enrichments/snpcontrolsets.rdata")
m<-match(bim$snp,f.all$SNP)
f.all<-f.all[m,]
df<-data.frame(bim_dt,f.all$closest450kdistance)

path="/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/maftssd/"
for (i in 1:23) {
cat(i,"\n")
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
df.chr<-df[which(df$chr==p),]
m<-match(r$V1,df.chr$pos)
r2<-data.frame(r,df.chr[m,c("GC_freq","CpG_freq","f.all.closest450kdistance")])
w<-which(is.na(r2$GC_freq))
write.table(r2,paste(path,"chr",i,".tmp",sep=""),sep=" ",quote=F,row.names=F,col.names=F)
}

pdf("GC_CpGcontent.pdf",height=6,width=6)
hist(bim_dt$GC_freq)
hist(bim_dt$CpG_freq)
dev.off()

#for i in `seq 1 23`; do
#echo $i
#mv chr$i chr$i.orig
#mv chr$i.tmp chr$i
#done


#bim_dt[,isGoDMC:=ifelse(cpgID%in%data$cpgname,TRUE,FALSE),]

#plot CG and CpG frequency for GoDMC cpgs and background 
#pdf("compare_seq_properties.pdf",height=3,width=4)
#ggplot(Illumina450_dt,aes(x=GC_freq,col=isGoDMC))+geom_density()
#ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
#ggplot(Illumina450_dt,aes(x=CpG_freq,col=isGoDMC))+geom_density()
#ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()
#dev.off()