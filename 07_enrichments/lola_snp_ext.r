message("lola snp ext")
library(LOLA)
library(dplyr)
library(GenomicRanges)
library(qvalue)

#load("../data/lola/snp_granges.rdata")

tfbsdb <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

message("global")

load("snpsetsforLOLA.RData")

mqtl_snps<-mqtlset
mqtl_snps$snpchr<-gsub("chr23","chrX",mqtl_snps$snpchr)

GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()

for (i in 1:10){
  controlset[[i]]$snpchr<-gsub("chr23","chrX",controlset[[i]]$snpchr)
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

#for (i in 1:10){
i=1
lola_res=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], tfbsdb, cores=5)
#plotLOLA(locResults_all=lola_res,plot_pref=paste0("snp_extbg_",i),height=6,width=9)

#}
save(lola_res, file="../results/lola_core_mqtlsnp.rdata")
####

regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")

anno<-regionDB$regionAnno
files<-anno[anno$collection%in%c("roadmap_epigenomics"),anno$filename]
#s<-seq(0,length(files),100)
#for (f in 1:length(s)){ 
#cat(f,"\n")
#getRegionSet(regionDB, collections="roadmap_epigenomics",filenames=files[(s[f]+1):s[f+1]])
getRegionSet(regionDB, collections="roadmap_epigenomics",filenames=files)
load("snpsetsforLOLA.RData")

mqtl_snps<-mqtlset
mqtl_snps$snpchr<-gsub("chr23","chrX",mqtl_snps$snpchr)

GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))

#use sampled backgrounds
lola_res_list=list()
snp_bg_gr_list=list()

for (i in 1:10){
  controlset[[i]]$snpchr<-gsub("chr23","chrX",controlset[[i]]$snpchr)
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

#for (i in 1:10){
i=1
lola_res=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], regionDB, cores=5)
#plotLOLA(locResults_all=lola_res,plot_pref=paste0("snp_extbg_",i),height=6,width=9)
#lola_res_list[[s]]<-lola_res
#}
#}
save(lola_res, file="../results/lola_ext_mqtlsnp.rdata")

#Data.table with enrichment results. Rows correspond to individual
#     pairwise fisher's tests comparing a single userSet with a single
#     databaseSet. The columns in this data.table are: userSet and
#     dbSet: index into their respective input region sets. pvalueLog:
#     -log10(pvalue) from the fisher's exact result; oddsRatio: result
#     from the fisher's exact test; support: number of regions in
#     userSet overlapping databaseSet; rnkPV, rnkLO, rnkSup: rank in
#     this table of p-value, oddsRatio, and Support respectively. The
#     -value is the negative natural log of the p-value returned from a
#     one-sided fisher's exact test. maxRnk, meanRnk: max and mean of
#     the 3 previous ranks, providing a combined ranking system. b, c,
#     d: 3 other values completing the 2x2 contingency table (with
#     support). The remaining columns describe the dbSet for the row.

#     If you have the qvalue package installed from bioconductor,
#     runLOLA will add a q-value transformation to provide FDR scores
#     automatically.
