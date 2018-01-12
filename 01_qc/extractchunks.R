path="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis"

load(paste0(path,"/results/16/16_clumped.rdata"))
nrow(clumped)
clumped<-clumped[clumped$pval<1e-14,]
nrow(clumped)
snps<-read.table(paste0(path,"/01_qc/indep.clump.prune.in"))

clumped<-clumped[which(clumped$snp%in%snps[,1]),]
nrow(clumped)

#a<-clumped
clumped$id<-as.character(paste(clumped$snp,clumped$cpg,sep="_"))
#a14.cis.out<-a[which(a$cis==F & a$pval<1e-14),]
#dim(a14.cis.out)
#[1] 156857     22 cis
# [1] 21526    29 trans
#[1] 23103    29
a14.cis.out<-clumped

a14.cis.out<-a14.cis.out[which(a14.cis.out$snpchr=="chr20"),] #332

a14.cis.out.filtered<-data.frame()

for (i in 1:length(unique(a14.cis.out$cpgchr))){
a14.cis.chr<-a14.cis.out[which(a14.cis.out$cpgchr==unique(a14.cis.out$cpgchr)[i]),]
o<-order(a14.cis.chr$cpgpos)
a14.cis.chr<-a14.cis.chr[o,]

cpgdiff<-NULL
for (j in 1:nrow(a14.cis.chr)-1){
cpgdiff[j]<-a14.cis.chr$cpgpos[j+1]-a14.cis.chr$cpgpos[j]
}

w<-which(cpgdiff<250000)
a14.cis.chr<-a14.cis.chr[-w,]
a14.cis.out.filtered<-rbind(a14.cis.out.filtered,a14.cis.chr)
}
#83

#o<-order(a14.cis.out$pval)
#a14.cis.out<-a14.cis.out[o,]
#probe<-unique(a14.cis.out$cpg)
#m<-match(probe,a14.cis.out$cpg)
#a14.cis.out<-a14.cis.out[m,]

#a14.cis.out<-a[which(a$cis==T & a$pval<1e-14&a$cpgchr=="chr20"),]
#dim(a14.cis.out)
#[1] 3213   21
#dim(data.frame(table(a14.cis.out$snp)))
#[1] 2497    2

#o<-order(as.numeric(as.character(a14.cis.out$cpgpos)))
#a14.cis.out<-a14.cis.out[o,]

chunks<-unique(data.frame(a14.cis.out$chunk))
chunks<-as.character(chunks[,1])
write.table(chunks,"chr20_chunks.txt",sep=" ",quote=F,row.names=F,col.names=F)

write.table(a14.cis.out$id,"chr20_ids.txt",sep=" ",quote=F,row.names=F,col.names=F)
