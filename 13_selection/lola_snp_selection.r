
message("lola snp ext")
library(data.table)
library(LOLA)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(qvalue)
library(gridExtra)
#load("../data/lola/snp_granges.rdata")
##load cell type conversion and colors
cellType_conversions=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/CellTypes.tsv",drop="collection")
colors=fread("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/color.tsv")

seldb <- loadRegionDB("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/regionDB2/hg19")

message("selection")

load("snpsetsforLOLA6prop.RData")


mqtl_snps<-mqtlset
GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()
for (i in 1:10){
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

for (i in 1:10){
lola_res=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], seldb, cores=5)
#plotLOLA(locResults_all=lola_res,plot_pref=paste0("snp_extbg_",i),height=6,width=9)
}

f.all.m<-rbind(mqtl_snps,controlset[[1]])
load("../results/enrichments/snpcontrolsetsGC_CpGcontent.rdata")
f.all<-r.all
m<-match(f.all.m$SNP,f.all$SNP)
f.all<-f.all[m,]
table(f.all$mQTL)

p1 <- ggplot(f.all, aes(MAF, fill=mQTL)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all, aes(nproxies, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all, aes(tssdist, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all,aes(closest450kdistance,fill=mQTL)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all, aes(CpG_freq, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all, aes(GC_freq, fill=mQTL)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensitymatched6prop.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()


load("snpsetsforLOLAtrans6prop.RData")
mqtl_snps<-mqtlset
GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()
for (i in 1:10){
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

for (i in 1:10){
lola_res_trans=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], seldb, cores=5)
}

load("snpsetsforLOLAcis6prop.RData")
mqtl_snps<-mqtlset
GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()
for (i in 1:10){
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

for (i in 1:10){
lola_res_cis=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], seldb, cores=5)
}


########DOESN'T WORK#######################################

seldb <- loadRegionDB("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/regionDB/hg19")

message("selection")

load("snpsetsforLOLA.RData")


mqtl_snps<-mqtlset
GoDMC_snp_gr=unique(with(mqtl_snps,GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
#use sampled backgrounds

snp_bg_gr_list=list()
for (i in 1:10){
  snp_bg_gr=unique(with(controlset[[i]],GRanges(seqnames = Rle(snpchr), IRanges(start=min, end=max),strand=Rle("*"))))
  snp_bg_gr_list[i]=unique(c(snp_bg_gr,GoDMC_snp_gr))
}

for (i in 1:10){
lola_res=runLOLA(GoDMC_snp_gr, snp_bg_gr_list[[i]], seldb, cores=5)
}


#testdata
library("LOLA")
regionDB <- loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

regionSetA = readBed("lola_vignette_data/setA_100.bed")
regionSetB = readBed("lola_vignette_data/setB_100.bed")
regionSetC = readBed("lola_vignette_data/setC_100.bed")
activeDHS = readBed("lola_vignette_data/activeDHS_universe.bed")
