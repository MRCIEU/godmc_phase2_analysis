# Library
library(GenomicRanges)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationHub)
library(data.table)
library(tidyverse)


arguments <- commandArgs(T)
i <- as.character(arguments[1])
var1 <- as.character(arguments[2])
var2 <- as.character(arguments[3])


print(i)
print(var1)
print(var2)


#i <-"chr18_chrX"
#var1 <- "chr18"
#var2 <- "chrX"



#1. Load and format data
#2. Select GRanges for SNP+proxies
#3. Select GRanges for CpG
#4. Select Granges for Hi-C data
#5. Overlap functions

#1.
rao <- read.table(paste0("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/",i,"/MAPQGE30/",i,".NORM"), header=F)
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
load("/panfs/panasas01/shared-godmc/1kg_reference_ph3/snpcontrolsets_selection.rdata") #284819, 24443 trans, 260376 cis
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE )) #271724, 23117 trans, 248607 cis

colnames(rao) <- c("bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe")
rao$interaction <- 1:nrow(rao)
rao$chr.oe<-gsub("chrX","chr23",rao$chr.oe)
rao$chr.bait<-gsub("chrX","chr23",rao$chr.bait)

var1 <-gsub("chrX", "chr23", var1)
var2 <-gsub("chrX", "chr23", var2)
print(var1)
print(var2)

clumped$code <- paste(clumped$cpg, clumped$snp)
trans_mqtl <- subset(clumped, clumped$cpgchr != clumped$snpchr) # 18584 trans mQTLs (interchromosomal)

#how many specific chr mQTL pairs are there?
trans_mqtl2 <-subset(trans_mqtl, (cpgchr ==var1 | cpgchr ==var2) & (snpchr==var1 | snpchr ==var2))
table(trans_mqtl2$cpgchr)
write.table(trans_mqtl2, file=paste0("trans_mqtl_clumped_",i,".tsv"), sep="\t", quote=F, row.names=F)



#2. take proxies, subset trans_mqtls for non duplicate snps, join non duplicated snps with proxies keeping matches only, GRanges
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
names(ldinfo) <- c("snp", "snplow", "snphigh", "snpproxies")
#temp <- subset(trans_mqtl2, !duplicated(snp))
temp <- inner_join(trans_mqtl2, ldinfo, "snp") 
snps <- GRanges(seqnames=temp$snpchr, ranges=IRanges(temp$snplow, temp$snphigh), strand="*")
names(snps) <- temp$snp # 12252 unique SNPs
mcols(snps) <- temp

#3. 
#temp <- subset(trans_mqtl2, !duplicated(cpg))
cpgs <- GRanges(seqnames=trans_mqtl2$cpgchr, ranges=IRanges(trans_mqtl2$cpgpos, trans_mqtl2$cpgpos), strand="*")
names(cpgs) <- trans_mqtl2$cpg # There are 15759 unique CpGs
mcols(cpgs) <- trans_mqtl2

#save(snps, cpgs, file=paste0("trans_mqtl_granges_",i,".rdata"))

#4. 
hic_bait <- GRanges(seqnames=rao$chr.bait, ranges=IRanges(rao$bait.start, rao$bait.end), strand="*")
names(hic_bait) <- rao$interaction
mcols(hic_bait) <- rao

hic_oe <- GRanges(seqnames=rao$chr.oe, ranges=IRanges(rao$oe.start, rao$oe.end), strand="*")
names(hic_oe) <- rao$interaction
mcols(hic_oe) <- rao

# lift over
#hub <- AnnotationHub()
#chain <- query(hub, "hg38ToHg19")[[1]]
#hic_bait <- liftOver(hic_bait, chain)
#hic_bait <- unlist(hic_bait)
#hic_oe <- liftOver(hic_oe, chain)
#hic_oe <- unlist(hic_oe)

#save(hic_bait, hic_oe, file=paste0("hi_c_ranges_",i,".rdata"))

#5.

bait.overlap.cpg <- function() {
  overlaps <- findOverlaps(hic_bait, cpgs, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$CpG <- names(cpgs[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp <- function() {
  overlaps <- findOverlaps(hic_bait, snps, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$SNP <- names(snps[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg <- function() {
  overlaps <- findOverlaps(hic_oe, cpgs, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$CpG <- names(cpgs[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp <- function() {
  overlaps <- findOverlaps(hic_oe, snps, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$SNP <- names(snps[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

bait_hic_cpg <- bait.overlap.cpg()
bait_hic_snp <- bait.overlap.snp()
oe_hic_cpg <- oe.overlap.cpg()
oe_hic_snp <- oe.overlap.snp()


# snp in bait
snp_in_bait <- merge(oe_hic_cpg, subset(bait_hic_snp, select=c(interaction, SNP)), by="interaction")
snp_in_bait$code <- paste(snp_in_bait$CpG, snp_in_bait$SNP)
snp_in_bait <- snp_in_bait[snp_in_bait$code %in% clumped$code, ]

# cpg in bait
cpg_in_bait <- merge(oe_hic_snp, subset(bait_hic_cpg, select=c(interaction, CpG)), by="interaction")
cpg_in_bait$code <- paste(cpg_in_bait$CpG, cpg_in_bait$SNP)
cpg_in_bait <- cpg_in_bait[cpg_in_bait$code %in% clumped$code, ]

data <- list(snp_in_bait, cpg_in_bait)
print(data)
save(data, file=paste0("data2_",i,".Rdata"))

#write.table(bait_hic_cpg, file=paste0("bait_hic_cpg_",i,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(bait_hic_snp, file=paste0("bait_hic_snp_",i,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(oe_hic_cpg, file=paste0("oe_hic_cpg_",i,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(oe_hic_snp, file=paste0("oe_hic_snp_",i,".tsv"), sep="\t", quote=F, row.names=F)
