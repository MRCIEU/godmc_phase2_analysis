library(data.table)
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicFeatures")
library(meffil)

arguments<-commandArgs(T)
i<-as.numeric(arguments[1])


#FUNCTION from MATT -really fast!
genomic.nearest <- function(positions, intervals) {
     stopifnot(all(c("chromosome","start","end") %in% colnames(intervals)))
     stopifnot(all(c("chromosome","position") %in% colnames(positions)))
 
     events <- rbind(data.frame(chromosome=intervals$chromosome,
                                position=intervals$start,
                                type="start",
                                id=1:nrow(intervals)),
                     data.frame(chromosome=intervals$chromosome,
                                position=intervals$end,
                                type="end",
                                id=1:nrow(intervals)),
                     data.frame(chromosome=positions$chromosome,
                                position=positions$position,
                                type="position",
                                id=1:nrow(positions)))
     events <- events[order(events$chromosome, events$position, decreasing=F),]
 
     before <- (1:nrow(events))-1
     before[1] <- NA
     for (position.idx in which(events$type == "position" & events$type[before] == "position"))
         if (position.idx > 1)
             before[position.idx] <- before[position.idx-1] ## 1s                                                                                                                                               
 
     after <- (1:nrow(events))+1
     after[length(after)] <- NA
     for(position.idx in rev(which(events$type == "position" & events$type[after] == "position")))
         if (position.idx < nrow(events))
             after[position.idx] <- after[position.idx+1]
 
 
     dist.before <- events$position - events$position[before]
     dist.before[which(events$chromosome != events$chromosome[before])] <- NA
     dist.after <- events$position[after] - events$position
     dist.after[which(events$chromosome != events$chromosome[after])] <- NA
     events$nearest <- NA
     idx <- which(!is.na(dist.before) & (is.na(dist.after) | dist.before <= dist.after))
     events$nearest[idx] <- before[idx]
     idx <- which(!is.na(dist.after) & (is.na(dist.before) | dist.before > dist.after))
     events$nearest[idx] <- after[idx]
 
     events$nearest <- events$id[events$nearest]
     positions$nearest <- NA
     idx <- which(events$type == "position")
     positions$nearest[events$id[idx]] <- events$nearest[idx]
     positions$nearest
 }

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/maftssd/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
names(r)<-c("pos","maf","tss","nproxies")
r$strand<-"+"
if (i<23) {r$chr<-p}
if (i==23) {r$chr<-"chrX"}

r_dt=as.data.table(r)
r_dt[,snpstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
r_dt[,snpend_pre:=ifelse(strand=="-",pos+500,pos+501),]

w<-which(r_dt$snpstart_pre<0)

if(length(w)>0){r_dt$snpstart_pre[w]<-1
r_dt$snpend_pre[w]<-500
}



#collapse overlaps
gr_range = with(r_dt,GRanges(seqnames=chr,ranges=IRanges(snpstart_pre,snpend_pre)))
gr_snp = with(r_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_snp, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

#
r_dt[,snpstart:=start(gr_range[overlap_red$subjectHit])]
r_dt[,snpend:=end(gr_range[overlap_red$subjectHit])]
r_dt[,NsubjectHits:=overlap_red$NsubjectHits]

hg19_gr=with(r_dt, GRanges(seqnames = Rle(chr), IRanges(start=snpstart, end=snpend),strand=Rle(strand),ID=pos))

seq_hg19=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_gr)

#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~ 3000 )
#Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]

r_dt[,GC_freq:=letterFrequency(seq_hg19, "CG", as.prob=T),]
r_dt[,CpG_freq:=dinucleotideFrequency(seq_hg19, step=2, as.prob=T)[,"CG"],]
#

#collapse overlaps

#check that center of seqeunce is always CpG (should be only the non CG probes and those that got merged into another region ~ 3000 )
#Illumina450_dt[NsubjectHits==1&subseq(seq_Illumina450,start=500,end=501)!="CG"]


r<-data.frame(r_dt)
##
y<-meffil.get.features("450k")

if(i==23){
j<-"X"
cpgpos<-y[which(y$chromosome==paste0("chr",j)),]}

if(i<23){cpgpos<-y[which(y$chromosome==paste0("chr",i)),]}
cpgpos$position<-as.numeric(as.character(cpgpos$position))

r$chromosome<-r$chr
r$position<-r$pos
cpgpos$start<-cpgpos$position
cpgpos$end<-cpgpos$position

res<-genomic.nearest(positions=r,intervals=cpgpos)
res<-cpgpos[res,c("name","position")]
cpgdist<-abs(res$position-r$pos)

r2<-data.frame(r$pos,r$maf,r$tss,r$nproxies,r$GC_freq,r$CpG_freq,cpg450kdist=cpgdist)

write.table(r2,paste(path,"chr",i,".tmp",sep=""),sep=" ",quote=F,row.names=F,col.names=F)


#for i in `seq 1 23`; do
#echo $i
##mv chr$i chr$i.orig
#mv chr$i chr$i.tmp
#done

#for i in `seq 1 23`; do
#echo $i
##mv chr$i chr$i.orig
#awk '{print $1, $2,$3,$4,$5,$6}' <chr$i.tmp >chr$i
#done

#bim_dt[,isGoDMC:=ifelse(cpgID%in%data$cpgname,TRUE,FALSE),]

#plot CG and CpG frequency for GoDMC cpgs and background 
#pdf("compare_seq_properties.pdf",height=3,width=4)
#ggplot(Illumina450_dt,aes(x=GC_freq,col=isGoDMC))+geom_density()
#ggplot(Illumina450_dt,aes(x=GC_freq))+geom_density()
#ggplot(Illumina450_dt,aes(x=CpG_freq,col=isGoDMC))+geom_density()
#ggplot(Illumina450_dt,aes(x=CpG_freq))+geom_density()
#dev.off()

