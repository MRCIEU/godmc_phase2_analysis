# setwd("I:/medewerkers/Koen/godmc/pchic/")
library(GenomicRanges)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationHub)
library(data.table)
library(tidyverse)
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


# Any interaction IDs that are the same on bait_pchic_cpg and oe_pchic_snp etc?

a <- merge(bait_pchic_cpg, subset(oe_pchic_snp, select=c(interaction, SNP)), by="interaction")
a$code <- paste(a$CpG, a$SNP)

b <- merge(oe_pchic_cpg, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction")
b$code <- paste(b$CpG, b$SNP)

clumped$code <- paste(clumped$cpg, clumped$snp)

a <- a[a$code %in% clumped$code, ]
b <- b[b$code %in% clumped$code, ]
ab <- rbind(a, b)

table(ab$oeChr == ab$baitChr)
table(a$SNP %in% clumped$snp)

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

# enrichment

reflection.fun <- function() {
  
  bait.overlap.cpg.reflection <- function(cpg) {
    cpg <- feats_reflection[names(feats_reflection) %in% cpg]
    overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
    result <- mcols(pchic_bait[queryHits(overlaps)])
    result$CpG <- names(cpg[subjectHits(overlaps)])
    return(result[order(result$CpG), ])
  }
  
  pchic_bait <- pchic_bait[na.omit(as.character(pchic_bait$baitChr) == as.character(pchic_bait$oeChr)), ]
  pchic_oe <- pchic_oe[na.omit(as.character(pchic_oe$baitChr) == as.character(pchic_oe$oeChr)), ]
  clumped <- clumped[na.omit(clumped$cpgchr == clumped$snpchr), ]
  
  # reflection
  clumped$reflection <- NA
  clumped$reflection[clumped$snppos < clumped$cpgpos] <- clumped$cpgpos[clumped$snppos < clumped$cpgpos] - clumped$snppos[clumped$snppos < clumped$cpgpos]
  clumped$reflection[clumped$snppos > clumped$cpgpos] <- clumped$cpgpos[clumped$snppos > clumped$cpgpos] + clumped$snppos[clumped$snppos > clumped$cpgpos]
  clumped <- clumped[!is.na(clumped$reflection), ]
  
  # circularize
  clumped <- lapply(unique(clumped$cpgchr), function(x) {
    clumped <- clumped[clumped$cpgchr == x, ]
    clumped$maxpos <- max(clumped$cpgpos, na.rm=T)
    clumped
  })
  clumped <- do.call(rbind, clumped)
  clumped$reflection[clumped$reflection < 0]  <- clumped$maxpos[clumped$reflection < 0] + clumped$reflection[clumped$reflection < 0]
  clumped$reflection[clumped$reflection > clumped$maxpos]  <- clumped$reflection[clumped$reflection > clumped$maxpos] - clumped$maxpos[clumped$reflection > clumped$maxpos]
  
  # add reflection positions to feats
  feats_reflection <- feats[match(clumped$cpg, names(feats))]
  start(feats_reflection) <- 1
  end(feats_reflection) <- max(end(feats_reflection))
  start(feats_reflection) <- clumped$reflection
  end(feats_reflection) <- clumped$reflection + 1
  
  bait_pchic_cpg_reflection <- bait.overlap.cpg.reflection(unique(clumped$cpg)) # reflection in promoter
  bait_pchic_cpg <- bait.overlap.cpg(unique(clumped$cpg)) # cpg in promoter
  oe_pchic_snp <- oe.overlap.snp(unique(clumped$snp)) # snp in interacting region
  
  a <- merge(bait_pchic_cpg, subset(oe_pchic_snp, select=c(interaction, SNP)), by="interaction") # cpg
  a$code <- paste(a$CpG, a$SNP)
  a <- a[a$code %in% clumped$code, ]
  
  b <- merge(bait_pchic_cpg_reflection, subset(oe_pchic_snp, select=c(interaction, SNP)), by="interaction") # reflection of cpg
  b$code <- paste(b$CpG, b$SNP)
  b <- b[b$code %in% clumped$code, ]
  
  # perform fisher test
  dat <- data.frame(c(length(unique(a$SNP)), length(unique(clumped$snp))), c(length(unique(b$SNP)), length(unique(clumped$snp))))
  fisher <- unlist(fisher.test(dat))[c(4, 1, 2, 3)]
  list(table=dat, test=fisher)
}

reflection <- reflection.fun()
print(reflection$test)

permutation.fun <- function(chr, n=10000) {
  set.seed(1)
  
  # only for cpg snp pairs on same chromosome, since pchic interactions almost only on same chr
  # seperate permutations per chr (fair distribution)
  pchic_bait <- pchic_bait[na.omit(as.character(pchic_bait$baitChr) == as.character(pchic_bait$oeChr) & as.character(pchic_bait$baitChr) == chr), ]
  pchic_oe <- pchic_oe[na.omit(as.character(pchic_oe$baitChr) == as.character(pchic_oe$oeChr) & as.character(pchic_oe$baitChr) == chr), ]
  clumped <- clumped[na.omit(clumped$cpgchr == clumped$snpchr & clumped$cpgchr == paste0("chr", chr)), ]

  bait_pchic_cpg <- bait.overlap.cpg(unique(clumped$cpg)) # cpg in promoter
  oe_pchic_snp <- oe.overlap.snp(unique(clumped$snp)) # snp in interacting region
  a <- merge(bait_pchic_cpg, subset(oe_pchic_snp, select=c(interaction, SNP)), by="interaction")
  a$code <- paste(a$CpG, a$SNP)
  len <- length(unique(a$SNP[a$code %in% clumped$code]))

  null <- c()
  for (i in 1:n) {
    # randomize interactions as null distribution
    clumped$cpg <- sample(clumped$cpg)
    clumped$code <- paste(clumped$cpg, clumped$snp)
    null <- c(null, length(unique(a$SNP[a$code %in% clumped$code])))
  }
  es <- mean((len + 1) / (null + 1))
  pval <- (sum(len <= null) + 1) / (length(null) + 1)
  data.frame(es, pval)
}

library(BiocParallel)
chrs <- unique(as.character(pchic_bait$baitChr))
chrs <- chrs[paste0("chr", chrs) %in% clumped$cpgchr]
permutation <- bplapply(chrs, permutation.fun, BPPARAM=MulticoreParam(length(chrs)))
enrichment <- do.call(rbind, permutation)
enrichment <- list(table=enrichment, test=apply(enrichment, 2, mean))
print(enrichment$test)

save(ab, abc, enrichment, reflection, file="../results/enrichments/mqtl_pchic.rdata")

