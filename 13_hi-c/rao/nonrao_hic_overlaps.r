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
pnum <- (arguments[4])

print(i)
print(var1)
print(var2)
print(pnum)

#i <-"chr15_chrX"
#var1 <- "chr15"
#var2 <- "chrX"
#pnum <-4


#1. Load and format data
#2. Select GRanges for SNP+proxies
#3. Select GRanges for cpg_samp
#4. Select Granges for Hi-C data
#5. Overlap functions
##6. Permutations (next script)
##7. Plots and tables (next script)

#1.
rao <- read.table(paste0("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/",i,"/MAPQGE30/",i,".NORM"), header=F)
load(paste0("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/permutations/permutations",pnum,".rdata"))
load("/panfs/panasas01/shared-godmc/1kg_reference_ph3/snpcontrolsets_selection.rdata") #284819, 24443 trans, 260376 cis

colnames(rao) <- c("bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe")
rao$interaction <- 1:nrow(rao)
rao$chr.oe<-gsub("chrX","chr23",rao$chr.oe)
rao$chr.bait<-gsub("chrX","chr23",rao$chr.bait)

var1 <-gsub("chrX", "chr23", var1)
var2 <-gsub("chrX", "chr23", var2)
print(var1)
print(var2)

trans_mqtl <- subset(perm2, perm2$cpgchr_samp != perm2$snpchr) # 18584 trans mQTLs (interchromosomal)

#how many specific chr mQTL pairs are there?
trans_mqtl2 <-subset(trans_mqtl, (cpgchr_samp ==var1 | cpgchr_samp ==var2) & (snpchr==var1 | snpchr ==var2))
table(trans_mqtl2$cpgchr_samp)
#write.table(trans_mqtl2, file=paste0("nontrans_mqtl_perm_",pnum,"_",i,".tsv"), sep="\t", quote=F, row.names=F)


#2. take proxies, subset trans_mqtls for non duplicate snps, join non duplicated snps with proxies keeping matches only, GRanges
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
names(ldinfo) <- c("snp", "snplow", "snphigh", "snpproxies")
temp <- subset(trans_mqtl2, !duplicated(snp))
temp <- inner_join(temp, ldinfo, "snp") 
snps <- GRanges(seqnames=temp$snpchr, ranges=IRanges(temp$snplow, temp$snphigh), strand="*")
names(snps) <- temp$snp # 12252 unique SNPs
mcols(snps) <- temp

#3. 
temp <- subset(trans_mqtl2, !duplicated(cpg_samp))
cpg_samps <- GRanges(seqnames=temp$cpgchr_samp, ranges=IRanges(temp$cpg_pos_samp, temp$cpg_pos_samp), strand="*")
names(cpg_samps) <- temp$cpg_samp # There are 15759 unique cpg_samps
mcols(cpg_samps) <- temp

#save(snps, cpg_samps, file=paste0("trans_mqtl_granges_",i,".rdata"))

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

#4.

bait.overlap.cpg_samp <- function() {
  overlaps <- findOverlaps(hic_bait, cpg_samps, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$cpg_samp <- names(cpg_samps[subjectHits(overlaps)])
  return(result[order(result$cpg_samp), ])
}

bait.overlap.snp <- function() {
  overlaps <- findOverlaps(hic_bait, snps, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$SNP <- names(snps[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg_samp <- function() {
  overlaps <- findOverlaps(hic_oe, cpg_samps, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$cpg_samp <- names(cpg_samps[subjectHits(overlaps)])
  return(result[order(result$cpg_samp), ])
}

oe.overlap.snp <- function() {
  overlaps <- findOverlaps(hic_oe, snps, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$SNP <- names(snps[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

bait_hic_cpg_samp <- bait.overlap.cpg_samp() # cpg_samp in bait 13556
bait_hic_snp <- bait.overlap.snp() # snp in bait 882383
oe_hic_cpg_samp <- oe.overlap.cpg_samp() # cpg_samp in interacting region 16388
oe_hic_snp <- oe.overlap.snp() # snp in interacting region 575885

# as data table as merge on null row (no overlap) datasets breaks the permutation loops

bait_hic_cpg_samp <- as.data.table(bait_hic_cpg_samp)
bait_hic_snp <- as.data.table(bait_hic_snp)
oe_hic_cpg_samp <- as.data.table(oe_hic_cpg_samp)
oe_hic_snp <- as.data.table(oe_hic_snp)


# snp in bait
snp_in_bait <- merge(oe_hic_cpg_samp, subset(bait_hic_snp, select=c(interaction, SNP)), by="interaction") # snp in bait
snp_in_bait$code <- paste(snp_in_bait$cpg_samp, snp_in_bait$SNP)
snp_in_bait <- snp_in_bait[snp_in_bait$code %in% perm2$code, ]
# cpg_samp in bait
cpg_in_bait <- merge(oe_hic_snp, subset(bait_hic_cpg_samp, select=c(interaction, cpg_samp)), by="interaction") # cpg_samp in bait
cpg_in_bait$code <- paste(cpg_in_bait$cpg_samp, cpg_in_bait$SNP)
cpg_in_bait <- cpg_in_bait[cpg_in_bait$code %in% perm2$code, ]

data <- list(snp_in_bait, cpg_in_bait)	
print(data)
save(data, file=paste0("nondata_",i,"_perm_",pnum,".Rdata"))

#write.table(bait_hic_cpg_samp, file=paste0("nonbait_hic_cpg_perm_",i,"_",pnum,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(bait_hic_snp, file=paste0("nonbait_hic_snp_perm_",i,"_",pnum,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(oe_hic_cpg_samp, file=paste0("nonoe_hic_cpg_perm_",i,"_",pnum,".tsv"), sep="\t", quote=F, row.names=F)
#write.table(oe_hic_snp, file=paste0("nonoe_hic_snp_perm_",i,"_",pnum,".tsv"), sep="\t", quote=F, row.names=F)

#to do
# check the results
# combine the overlaps
# permutations using broken links and 1000 tests
# enrichment of overlaps in LOLA: background set to be non associated mQTLS that are in contact points. 
# enrichment of SNPs in GARFIELD: will take a bit of time. 
# update document for meeting on Wednesday


