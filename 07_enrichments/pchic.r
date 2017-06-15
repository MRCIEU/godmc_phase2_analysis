# setwd("I:/medewerkers/Koen/godmc/pchic/")
library(GenomicRanges)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationHub)
library(data.table)
feats <- features(FDb.InfiniumMethylation.hg19)
pchic <- fread("zcat ../data/misc/PCHiC_peak_matrix_cutoff5.tsv.gz") #http://dx.doi.org/10.1016/j.cell.2016.09.037
pchic$interaction <- 1:nrow(pchic)

load("../results/16/16_clumped.rdata")

pchic_bait <- GRanges(seqnames=paste0("chr", pchic$baitChr), ranges=IRanges(start=pchic$baitStart, end=pchic$baitEnd, names=pchic$interaction))
mcols(pchic_bait) <- pchic
pchic_oe <- GRanges(seqnames=paste0("chr", pchic$oeChr), ranges=IRanges(start=pchic$oeStart, end=pchic$oeEnd, names=pchic$interaction))
mcols(pchic_oe) <- pchic

hub <- AnnotationHub()
chain <- query(hub, "hg38ToHg19")[[1]]

pchic_bait <- liftOver(pchic_bait, chain)
pchic_bait <- unlist(pchic_bait)

pchic_oe <- liftOver(pchic_oe, chain)
pchic_oe <- unlist(pchic_oe)

bait.overlap.cpg <- function(cpg) {
  cpg <- feats[names(feats) %in% cpg]
  overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

bait.overlap.snp <- function(snp) {
  split <- do.call(rbind, t(strsplit(snp, split=":")))
  chr <- split[, 1]
  loc <- as.numeric(split[, 2])
  snp <- GRanges(seqnames=chr, ranges=IRanges(start=loc, end=loc, names=snp))
  overlaps <- findOverlaps(pchic_bait, snp, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

oe.overlap.cpg <- function(cpg) {
  cpg <- feats[names(feats) %in% cpg]
  overlaps <- findOverlaps(pchic_oe, cpg, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

oe.overlap.snp <- function(snp) {
  split <- do.call(rbind, t(strsplit(snp, split=":")))
  chr <- split[, 1]
  loc <- as.numeric(split[, 2])
  snp <- GRanges(seqnames=chr, ranges=IRanges(start=loc, end=loc, names=snp))
  overlaps <- findOverlaps(pchic_oe, snp, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

bait_pchic_cpg <- bait.overlap.cpg(unique(clumped$cpg)) # cpg in promoter
bait_pchic_snp <- bait.overlap.snp(unique(clumped$snp)) # snp in promoter

oe_pchic_cpg <- oe.overlap.cpg(unique(clumped$cpg)) # cpg in interacting region
oe_pchic_snp <- oe.overlap.snp(unique(clumped$snp)) # snp in interacting region


# Any interaction IDs that are the same on bait_pchic_cpg and oe_pchic_snp?
library(dplyr)

a <- merge(bait_pchic_cpg, subset(oe_pchic_snp, select=c(interaction, SNP)), by="interaction")
a$code <- paste(a$CpG, a$SNP)

b <- merge(oe_pchic_cpg, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction")
b$code <- paste(b$CpG, b$SNP)

clumped$code <- paste(clumped$cpg, clumped$snp)

a <- a[a$code %in% clumped$code, ]
b <- b[b$code %in% clumped$code, ]
ab <- rbind(a, b)

table(ab$oeChr == ab$baitChr)
table(a$SNP %in% clumped$SNP)

abl <- subset(ab, select=c(code, Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, tB)) %>%
  as.data.frame %>%
  gather(key="celltype", value="chicagoscore", -code) %>%
  filter(chicagoscore > 5)

abc <- group_by(abl, code) %>%
  summarise(celltype=paste(celltype, collapse=","))

group_by(abl, celltype) %>%
  summarise(s=sum(chicagoscore))

# write.table(bait_pchic_cpg, file="bait_pchic_cpg.tsv", sep="\t", quote=F, row.names=F)
# write.table(bait_pchic_snp, file="bait_pchic_snp.tsv", sep="\t", quote=F, row.names=F)

# write.table(oe_pchic_cpg, file="oe_pchic_cpg.tsv", sep="\t", quote=F, row.names=F)
# write.table(oe_pchic_snp, file="oe_pchic_snp.tsv", sep="\t", quote=F, row.names=F)

save(ab, abc, file="../results/enrichments/mqtl_pchic.rdata")
