# Grubert Hi-C Data Analysis

# Library
library(GenomicRanges)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationHub)
library(data.table)
library(tidyverse)

# Load Data
load("../data/misc/Grubert_HiC_0.4_clean.Rdata")
load("../results/16/16_clumped.rdata")


# Format Data
hic$interaction <- 1:nrow(hic)
clumped$code <- paste(clumped$cpg, clumped$snp)

hic_bait <- GRanges(seqnames=paste0("chr", hic$chr), ranges=IRanges(start=hic$i.start, end=hic$i.end, names=hic$interaction))
mcols(hic_bait) <- hic
hic_oe <- GRanges(seqnames=paste0("chr", hic$chr), ranges=IRanges(start=hic$j.start, end=hic$j.end, names=hic$interaction))
mcols(hic_oe) <- hic


bait.overlap.cpg <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(hic_bait, cpg, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(hic_bait, snp, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(hic_oe, cpg, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(hic_oe, snp, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

bait.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(hic_bait, cpg, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp.reflection <- function(snp) {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(hic_bait, snp, ignore.strand=T)
  result <- mcols(hic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(hic_oe, cpg, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp.reflection <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(hic_oe, snp, ignore.strand=T)
  result <- mcols(hic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# cis
clumped <- subset(clumped, cpgchr == snpchr & abs(cpgpos - snppos) <= 1000000)
hic_oe <- hic_oe[na.omit(as.character(hic_oe$chr) == as.character(hic_oe$chr) & min(abs(c(hic_oe$j.end - hic_oe$i.start, hic_oe$j.start - hic_oe$i.end))) <= 1000000), ]
hic_bait <- hic_bait[na.omit(as.character(hic_bait$chr) == as.character(hic_bait$chr) & min(abs(c(hic_bait$j.end - hic_bait$i.start, hic_bait$j.start - hic_bait$i.end))) <= 1000000), ]

# cpgreflection
clumped$cpgreflection <- NA
clumped$cpgreflection[clumped$snppos < clumped$cpgpos] <- clumped$cpgpos[clumped$snppos < clumped$cpgpos] - clumped$snppos[clumped$snppos < clumped$cpgpos]
clumped$cpgreflection[clumped$snppos > clumped$cpgpos] <- clumped$cpgpos[clumped$snppos > clumped$cpgpos] + clumped$snppos[clumped$snppos > clumped$cpgpos]

# snpreflection
clumped$snpreflection <- NA
clumped$snpreflection[clumped$cpgpos < clumped$snppos] <- clumped$snppos[clumped$cpgpos < clumped$snppos] - clumped$cpgpos[clumped$cpgpos < clumped$snppos]
clumped$snpreflection[clumped$cpgpos > clumped$snppos] <- clumped$snppos[clumped$cpgpos > clumped$snppos] + clumped$cpgpos[clumped$cpgpos > clumped$snppos]

# outside chrom
clumped <- lapply(unique(clumped$cpgchr), function(x) {
  clumped <- clumped[clumped$cpgchr == x, ]
  clumped$cpgmax <- max(clumped$cpgpos, na.rm=T)
  clumped$snpmax <- max(clumped$snppos, na.rm=T)
  clumped
})
clumped <- do.call(rbind, clumped)

clumped$cpgreflection[clumped$cpgreflection < 0]  <- NA
clumped$cpgreflection[clumped$cpgreflection > clumped$cpgmax]  <- NA

clumped$snpreflection[clumped$snpreflection < 0]  <- NA
clumped$snpreflection[clumped$snpreflection > clumped$snpmax]  <- NA

clumped_cpg <- clumped[!is.na(clumped$cpgreflection), ]
clumped_snp <- clumped[!is.na(clumped$snpreflection), ]

# overlaps
bait_hic_cpg <- bait.overlap.cpg() # cpg in bait
bait_hic_snp <- bait.overlap.snp() # snp in bait
oe_hic_cpg <- oe.overlap.cpg() # cpg in interacting region
oe_hic_snp <- oe.overlap.snp() # snp in interacting region
bait_hic_cpg_reflection <- bait.overlap.cpg.reflection() # cpg reflection in bait
bait_hic_snp_reflection <- bait.overlap.snp.reflection() # snp reflection in bait
oe_hic_cpg_reflection <- oe.overlap.cpg.reflection() # cpg reflection in interacting region
oe_hic_snp_reflection <- oe.overlap.snp.reflection() # snp reflection in interacting region

# snp in bait
snp_in_bait <- merge(oe_hic_cpg, subset(bait_hic_snp, select=c(interaction, SNP)), by="interaction") # snp in bait
snp_in_bait$code <- paste(snp_in_bait$CpG, snp_in_bait$SNP)
snp_in_bait <- snp_in_bait[snp_in_bait$code %in% clumped$code, ]

snp_in_bait_reflection <- merge(oe_hic_cpg_reflection, subset(bait_hic_snp, select=c(interaction, SNP)), by="interaction") # snp in bait reflection
snp_in_bait_reflection$code <- paste(snp_in_bait_reflection$CpG, snp_in_bait_reflection$SNP)
snp_in_bait_reflection <- snp_in_bait_reflection[snp_in_bait_reflection$code %in% clumped$code, ]

# cpg in bait
cpg_in_bait <- merge(oe_hic_snp, subset(bait_hic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in bait
cpg_in_bait$code <- paste(cpg_in_bait$CpG, cpg_in_bait$SNP)
cpg_in_bait <- cpg_in_bait[cpg_in_bait$code %in% clumped$code, ]

cpg_in_bait_reflection <- merge(oe_hic_snp_reflection, subset(bait_hic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in bait reflection
cpg_in_bait_reflection$code <- paste(cpg_in_bait_reflection$CpG, cpg_in_bait_reflection$SNP)
cpg_in_bait_reflection <- cpg_in_bait_reflection[cpg_in_bait_reflection$code %in% clumped$code, ]

data <- list(snp_in_bait, cpg_in_bait)

library(ggplot2)

plotdata <- data.frame(data=c(rep("clumped", nrow(clumped)), rep("hic", nrow(hic)), rep("snp_in_bait", nrow(snp_in_bait)), rep("cpg_in_bait", nrow(cpg_in_bait))), distance=c(clumped$cpgpos - clumped$snppos, hic$dist, snp_in_bait$dist, cpg_in_bait$dist))
pdf("../results/enrichments/densityplot_grubert.pdf")
ggplot(plotdata, aes(distance, fill=data, alpha=0.5)) + geom_density(color=NA) + guides(alpha=F) + xlim(-1000000, 1000000) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.key = element_blank(), strip.background = element_rect(colour=NA, fill=NA), panel.border=element_blank(), panel.spacing = unit(2, "lines"))
dev.off()

# perform fisher test
dat <- data.frame(c(length(snp_in_bait$SNP), length(bait_hic_snp$SNP[!bait_hic_snp$SNP %in% snp_in_bait$SNP])), c(length(snp_in_bait_reflection$SNP), length(bait_hic_snp$SNP[!bait_hic_snp$SNP %in% snp_in_bait_reflection$SNP])))
fisher <- unlist(fisher.test(dat))[c(4, 1, 2, 3)]
snp_in_bait <- list(table=dat, test=fisher, unique_snps=length(unique(snp_in_bait$SNP)), percentage=length(unique(snp_in_bait$SNP)) / length(unique(bait_hic_snp$SNP)[!unique(bait_hic_snp$SNP) %in% unique(snp_in_bait$SNP)]))

dat <- data.frame(c(length(cpg_in_bait$CpG), length(bait_hic_cpg$CpG[!bait_hic_cpg$CpG %in% cpg_in_bait$CpG])), c(length(cpg_in_bait_reflection$CpG), length(bait_hic_cpg$CpG[!bait_hic_cpg$CpG %in% cpg_in_bait_reflection$CpG])))
fisher <- unlist(fisher.test(dat))[c(4, 1, 2, 3)]
cpg_in_bait <- list(table=dat, test=fisher, unique_cpgs=length(unique(cpg_in_bait$CpG)), percentage=length(unique(cpg_in_bait$CpG)) / length(unique(bait_hic_cpg$CpG)[!unique(bait_hic_cpg$CpG) %in% unique(cpg_in_bait$CpG)]))

enrichment <- list(snp_in_bait=snp_in_bait, cpg_in_bait=cpg_in_bait)
print(enrichment)

save(data, enrichment, file= "../results/enrichments/mqtl_hic_grubert.rdata")