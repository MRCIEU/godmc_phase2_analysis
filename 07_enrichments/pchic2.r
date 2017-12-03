# setwd("I:/medewerkers/Koen/godmc/pchic/")
library(GenomicRanges)
library(AnnotationHub)
library(FDb.InfiniumMethylation.hg19)
library(data.table)
library(tidyverse)
pchic <- fread("zcat ../data/misc/PCHiC_peak_matrix_cutoff5.tsv.gz") #http://dx.doi.org/10.1016/j.cell.2016.09.037
load("../results/16/16_clumped.rdata")

pchic$interaction <- 1:nrow(pchic)
clumped$code <- paste(clumped$cpg, clumped$snp)

pchic_bait <- GRanges(seqnames=paste0("chr", pchic$baitChr), ranges=IRanges(start=pchic$baitStart, end=pchic$baitEnd, names=pchic$interaction))
mcols(pchic_bait) <- pchic
pchic_oe <- GRanges(seqnames=paste0("chr", pchic$oeChr), ranges=IRanges(start=pchic$oeStart, end=pchic$oeEnd, names=pchic$interaction))
mcols(pchic_oe) <- pchic

bait.overlap.cpg <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(pchic_bait, snp, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_oe, cpg, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(pchic_oe, snp, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

bait.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp.reflection <- function(snp) {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(pchic_bait, snp, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_oe, cpg, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp.reflection <- function() {
  clumped <- clumped_snp
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(pchic_oe, snp, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# cis
clumped <- subset(clumped, cpgchr == snpchr & abs(cpgpos - snppos) <= 1000000)
pchic_oe <- pchic_oe[na.omit(as.character(pchic_oe$oeChr) == as.character(pchic_oe$baitChr) & min(abs(c(pchic_oe$oeEnd - pchic_oe$baitStart, pchic_oe$oeStart - pchic_oe$baitEnd))) <= 1000000), ]
pchic_bait <- pchic_bait[na.omit(as.character(pchic_bait$oeChr) == as.character(pchic_bait$baitChr) & min(abs(c(pchic_bait$oeEnd - pchic_bait$baitStart, pchic_bait$oeStart - pchic_bait$baitEnd))) <= 1000000), ]

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
bait_pchic_cpg <- bait.overlap.cpg() # cpg in promoter
bait_pchic_snp <- bait.overlap.snp() # snp in promoter
oe_pchic_cpg <- oe.overlap.cpg() # cpg in interacting region
oe_pchic_snp <- oe.overlap.snp() # snp in interacting region
bait_pchic_cpg_reflection <- bait.overlap.cpg.reflection() # cpg reflection in promoter
bait_pchic_snp_reflection <- bait.overlap.snp.reflection() # snp reflection in promoter
oe_pchic_cpg_reflection <- oe.overlap.cpg.reflection() # cpg reflection in interacting region
oe_pchic_snp_reflection <- oe.overlap.snp.reflection() # snp reflection in interacting region

# snp in promoter
snp_in_promoter <- merge(oe_pchic_cpg, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction") # snp in promoter
snp_in_promoter$code <- paste(snp_in_promoter$CpG, snp_in_promoter$SNP)
snp_in_promoter <- snp_in_promoter[snp_in_promoter$code %in% clumped$code, ]

snp_in_promoter_reflection <- merge(oe_pchic_cpg_reflection, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction") # snp in promoter reflection
snp_in_promoter_reflection$code <- paste(snp_in_promoter_reflection$CpG, snp_in_promoter_reflection$SNP)
snp_in_promoter_reflection <- snp_in_promoter_reflection[snp_in_promoter_reflection$code %in% clumped$code, ]

# cpg in promoter
cpg_in_promoter <- merge(oe_pchic_snp, subset(bait_pchic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in promoter
cpg_in_promoter$code <- paste(cpg_in_promoter$CpG, cpg_in_promoter$SNP)
cpg_in_promoter <- cpg_in_promoter[cpg_in_promoter$code %in% clumped$code, ]

cpg_in_promoter_reflection <- merge(oe_pchic_snp_reflection, subset(bait_pchic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in promoter reflection
cpg_in_promoter_reflection$code <- paste(cpg_in_promoter_reflection$CpG, cpg_in_promoter_reflection$SNP)
cpg_in_promoter_reflection <- cpg_in_promoter_reflection[cpg_in_promoter_reflection$code %in% clumped$code, ]

data <- list(snp_in_promoter, cpg_in_promoter)

library(ggplot2)

plotdata <- data.frame(data=c(rep("clumped", nrow(clumped)), rep("pchic", nrow(pchic)), rep("snp_in_promoter", nrow(snp_in_promoter)), rep("cpg_in_promoter", nrow(cpg_in_promoter))), distance=c(clumped$cpgpos - clumped$snppos, pchic$dist, snp_in_promoter$dist, cpg_in_promoter$dist))
pdf("../results/enrichments/densityplot.pdf")
ggplot(plotdata, aes(distance, fill=data, alpha=0.5)) + geom_density(color=NA) + guides(alpha=F) + xlim(-1000000, 1000000) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.key = element_blank(), strip.background = element_rect(colour=NA, fill=NA), panel.border=element_blank(), panel.spacing = unit(2, "lines"))
dev.off()

pdf("../results/enrichments/snp_in_promoter.pdf")
heatmap(as.matrix(snp_in_promoter[, 13:29]), labRow=NA, col=colorRampPalette(c("white", "red"))(256))
dev.off()

pdf("../results/enrichments/cpg_in_promoter.pdf")
heatmap(as.matrix(cpg_in_promoter[, 13:29]), labRow=NA, col=colorRampPalette(c("white", "red"))(256))
dev.off()

# perform fisher test
dat <- data.frame(c(length(snp_in_promoter$SNP), length(bait_pchic_snp$SNP[!bait_pchic_snp$SNP %in% snp_in_promoter$SNP])), c(length(snp_in_promoter_reflection$SNP), length(bait_pchic_snp$SNP[!bait_pchic_snp$SNP %in% snp_in_promoter_reflection$SNP])))
fisher <- unlist(fisher.test(dat))[c(4, 1, 2, 3)]
snp_in_promoter <- list(table=dat, test=fisher, unique_snps=length(unique(snp_in_promoter$SNP)), percentage=length(unique(snp_in_promoter$SNP)) / length(unique(bait_pchic_snp$SNP)[!unique(bait_pchic_snp$SNP) %in% unique(snp_in_promoter$SNP)]))

dat <- data.frame(c(length(cpg_in_promoter$CpG), length(bait_pchic_cpg$CpG[!bait_pchic_cpg$CpG %in% cpg_in_promoter$CpG])), c(length(cpg_in_promoter_reflection$CpG), length(bait_pchic_cpg$CpG[!bait_pchic_cpg$CpG %in% cpg_in_promoter_reflection$CpG])))
fisher <- unlist(fisher.test(dat))[c(4, 1, 2, 3)]
cpg_in_promoter <- list(table=dat, test=fisher, unique_cpgs=length(unique(cpg_in_promoter$CpG)), percentage=length(unique(cpg_in_promoter$CpG)) / length(unique(bait_pchic_cpg$CpG)[!unique(bait_pchic_cpg$CpG) %in% unique(cpg_in_promoter$CpG)]))

enrichment <- list(snp_in_promoter=snp_in_promoter, cpg_in_promoter=cpg_in_promoter)
print(enrichment)

save(data, enrichment, file="../results/enrichments/mqtl_pchic2.rdata")
