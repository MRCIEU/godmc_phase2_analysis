library(GenomicRanges)

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19/encode_tfbs/encode_tfbs.RData")
length(ret) #689
gr.list<-ret

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19/encode_tfbs/encode_tfbs_files.RData")
dim(ret)
#[1] 689

for (i in 1:length(gr.list)){
cat(i,"\n")
gr<-gr.list[[i]]
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  antibody=c(rep(ret[i,"antibody"], length(gr))),
  cellType=c(rep(ret[i,"cellType"], length(gr))),
  strands=strand(gr))

write.table(df, file=paste0("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs/tfbs_",ret$filename[i],".bed"), quote=F, sep="\t", row.names=F, col.names=F)
}
#0 based
# 50070
