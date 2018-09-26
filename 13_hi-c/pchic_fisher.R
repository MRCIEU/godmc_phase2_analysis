options(stringsAsFactors = F)
set.seed(1)

library(GenomicRanges)
library(BiocParallel)
library(data.table)

# load data
pchic <- fread("zcat ../data/misc/PCHiC_peak_matrix_cutoff5.tsv.gz") #http://dx.doi.org/10.1016/j.cell.2016.09.037
load("../results/16/16_clumped.rdata")
pchic <- as.data.frame(pchic)

# only cis-mqtls and only pchic on same chr
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
clumped <- clumped[clumped$cis == T, ]
pchic <- pchic[pchic$baitChr == pchic$oeChr, ]

# create unique interaction ID
pchic$interaction <- 1:nrow(pchic)

# create unique cpg-snp ID
clumped$code <- paste(clumped$cpg, clumped$snp)

# granges of promoters
pchic_bait <- GRanges(seqnames=paste0("chr", pchic$baitChr), ranges=IRanges(start=pchic$baitStart, end=pchic$baitEnd, names=pchic$interaction))
mcols(pchic_bait) <- pchic

# granges of interacting regions
pchic_oe <- GRanges(seqnames=paste0("chr", pchic$oeChr), ranges=IRanges(start=pchic$oeStart, end=pchic$oeEnd, names=pchic$interaction))
mcols(pchic_oe) <- pchic

# calculate overlap between promoter and cpg
bait.overlap.cpg <- function() {
  clumped <- clumped[!duplicated(clumped$cpg), ]
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

# calculate overlap between promoter and snp
bait.overlap.snp <- function() {
  clumped <- clumped[!duplicated(clumped$snp), ]
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(pchic_bait, snp, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# calculate overlap between interacting region and cpg
oe.overlap.cpg <- function() {
  clumped <- clumped[!duplicated(clumped$cpg), ]
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgpos, end=clumped$cpgpos, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_oe, cpg, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

# calculate overlap between interacting region and snp
oe.overlap.snp <- function() {
  clumped <- clumped[!duplicated(clumped$snp), ]
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snppos, end=clumped$snppos, names=clumped$snp))
  overlaps <- findOverlaps(pchic_oe, snp, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# calculate overlap between promoter and cpg reflection
bait.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg[!duplicated(clumped_cpg$cpg), ]
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_bait, cpg, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

# calculate overlap between promoter and snp reflection
bait.overlap.snp.reflection <- function(snp) {
  clumped <- clumped_snp[!duplicated(clumped_snp$snp), ]
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(pchic_bait, snp, ignore.strand=T)
  result <- mcols(pchic_bait[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# calculate overlap between interacting region and cpg reflection
oe.overlap.cpg.reflection <- function() {
  clumped <- clumped_cpg[!duplicated(clumped_cpg$cpg), ]
  cpg <- GRanges(seqnames=clumped$cpgchr, ranges=IRanges(start=clumped$cpgreflection, end=clumped$cpgreflection, names=clumped$cpg))
  overlaps <- findOverlaps(pchic_oe, cpg, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$CpG <- names(cpg[subjectHits(overlaps)])
  return(result[order(result$CpG), ])
}

# calculate overlap between interacting region and snp reflection
oe.overlap.snp.reflection <- function() {
  clumped <- clumped_snp[!duplicated(clumped_snp$snp), ]
  snp <- GRanges(seqnames=clumped$snpchr, ranges=IRanges(start=clumped$snpreflection, end=clumped$snpreflection, names=clumped$snp))
  overlaps <- findOverlaps(pchic_oe, snp, ignore.strand=T)
  result <- mcols(pchic_oe[queryHits(overlaps)])
  result$SNP <- names(snp[subjectHits(overlaps)])
  return(result[order(result$SNP), ])
}

# calculate cpg reflection
clumped$cpgreflection <- NA
clumped$cpgreflection[clumped$snppos < clumped$cpgpos] <- clumped$cpgpos[clumped$snppos < clumped$cpgpos] - clumped$snppos[clumped$snppos < clumped$cpgpos]
clumped$cpgreflection[clumped$snppos > clumped$cpgpos] <- clumped$cpgpos[clumped$snppos > clumped$cpgpos] + clumped$snppos[clumped$snppos > clumped$cpgpos]

# calculate snp reflection
clumped$snpreflection <- NA
clumped$snpreflection[clumped$cpgpos < clumped$snppos] <- clumped$snppos[clumped$cpgpos < clumped$snppos] - clumped$cpgpos[clumped$cpgpos < clumped$snppos]
clumped$snpreflection[clumped$cpgpos > clumped$snppos] <- clumped$snppos[clumped$cpgpos > clumped$snppos] + clumped$cpgpos[clumped$cpgpos > clumped$snppos]

# remove regions outside chromosome
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

# proportion in gene
genesymbol <- function(chr, loc) {
  
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  granges <- GRanges(seqnames = chr, ranges = IRanges(start = loc, end = loc))
  genes <- unlist(cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene"))
  near <- genes[nearest(granges, genes, ignore.strand = T)]
  select(org.Hs.eg.db, keys = names(near), columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
  
}

clumped_reflection <- clumped[!is.na(clumped$snpreflection) & !is.na(clumped$cpgreflection) & !clumped$snpchr == "chr23", ]
clumped_reflection$snpgenereflection <- genesymbol(clumped_reflection$snpchr, clumped_reflection$snpreflection)
clumped_reflection$cpggenereflection <- genesymbol(clumped_reflection$cpgchr, clumped_reflection$cpgreflection)
clumped_reflection$snpgene <- genesymbol(clumped_reflection$snpchr, clumped_reflection$snppos)
clumped_reflection$cpggene <- genesymbol(clumped_reflection$cpgchr, clumped_reflection$cpgpos)

print("Proportion in gene: mQTL snp and cpg")
print(prop.table(table(clumped_reflection$snpgene == clumped_reflection$cpggene)))
print("Proportion in gene: mQTL snp and cpg reflection")
print(prop.table(table(clumped_reflection$snpgene == clumped_reflection$cpggenereflection)))
print("Proportion in gene: mQTL snp reflection and cpg")
print(prop.table(table(clumped_reflection$snpgenereflection == clumped_reflection$cpggene))) 

# calculate overlaps
bait_pchic_cpg <- bait.overlap.cpg() # cpg in promoter
bait_pchic_snp <- bait.overlap.snp() # snp in promoter
oe_pchic_cpg <- oe.overlap.cpg() # cpg in interacting region
oe_pchic_snp <- oe.overlap.snp() # snp in interacting region
bait_pchic_cpg_reflection <- bait.overlap.cpg.reflection() # cpg reflection in promoter
bait_pchic_snp_reflection <- bait.overlap.snp.reflection() # snp reflection in promoter
oe_pchic_cpg_reflection <- oe.overlap.cpg.reflection() # cpg reflection in interacting region
oe_pchic_snp_reflection <- oe.overlap.snp.reflection() # snp reflection in interacting region

# snp in promoter & cpg in interacting region
snp_in_promoter <- merge(oe_pchic_cpg, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction") # snp in promoter
snp_in_promoter$code <- paste(snp_in_promoter$CpG, snp_in_promoter$SNP)
snp_in_promoter_code <- snp_in_promoter[snp_in_promoter$code %in% clumped$code, ] # is mqtl
snp_in_promoter_nocode <- snp_in_promoter[!snp_in_promoter$code %in% clumped$code, ] # is not mqtl

snp_in_promoter_reflection <- merge(oe_pchic_cpg_reflection, subset(bait_pchic_snp, select=c(interaction, SNP)), by="interaction") # snp in promoter reflection
snp_in_promoter_reflection$code <- paste(snp_in_promoter_reflection$CpG, snp_in_promoter_reflection$SNP)
snp_in_promoter_code_reflection <- snp_in_promoter_reflection[snp_in_promoter_reflection$code %in% clumped$code, ]
snp_in_promoter_nocode_reflection <- snp_in_promoter_reflection[!snp_in_promoter_reflection$code %in% clumped$code, ]

# cpg in promoter & snp in interacting region
cpg_in_promoter <- merge(oe_pchic_snp, subset(bait_pchic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in promoter
cpg_in_promoter$code <- paste(cpg_in_promoter$CpG, cpg_in_promoter$SNP)
cpg_in_promoter_code <- cpg_in_promoter[cpg_in_promoter$code %in% clumped$code, ]
cpg_in_promoter_nocode <- cpg_in_promoter[!cpg_in_promoter$code %in% clumped$code, ]

cpg_in_promoter_reflection <- merge(oe_pchic_snp_reflection, subset(bait_pchic_cpg, select=c(interaction, CpG)), by="interaction") # cpg in promoter reflection
cpg_in_promoter_reflection$code <- paste(cpg_in_promoter_reflection$CpG, cpg_in_promoter_reflection$SNP)
cpg_in_promoter_code_reflection <- cpg_in_promoter_reflection[cpg_in_promoter_reflection$code %in% clumped$code, ]
cpg_in_promoter_nocode_reflection <- cpg_in_promoter_reflection[!cpg_in_promoter_reflection$code %in% clumped$code, ]

# merge snp in promoter and cpg in promoter
snp_or_cpg_in_promoter_code <- rbind(snp_in_promoter_code, cpg_in_promoter_code[, c(1:31, 33, 32, 34)])
snp_or_cpg_in_promoter_code_reflection <- rbind(snp_in_promoter_code_reflection, cpg_in_promoter_code_reflection[, c(1:31, 33, 32, 34)])

snp_or_cpg_in_promoter_code <- snp_or_cpg_in_promoter_code[!duplicated(paste0(snp_or_cpg_in_promoter_code$code, snp_or_cpg_in_promoter_code$interaction)), ]
snp_or_cpg_in_promoter_code_reflection <- snp_or_cpg_in_promoter_code_reflection[!duplicated(paste0(snp_or_cpg_in_promoter_code_reflection$code, snp_or_cpg_in_promoter_code_reflection$interaction)), ]

snp_or_cpg_in_promoter_nocode <- rbind(snp_in_promoter_nocode, cpg_in_promoter_nocode[, c(1:31, 33, 32, 34)])
snp_or_cpg_in_promoter_nocode_reflection <- rbind(snp_in_promoter_nocode_reflection, cpg_in_promoter_nocode_reflection[, c(1:31, 33, 32, 34)])

snp_or_cpg_in_promoter_nocode <- snp_or_cpg_in_promoter_nocode[!duplicated(paste0(snp_or_cpg_in_promoter_nocode$code, snp_or_cpg_in_promoter_nocode$interaction)), ]
snp_or_cpg_in_promoter_nocode_reflection <- snp_or_cpg_in_promoter_nocode_reflection[!duplicated(paste0(snp_or_cpg_in_promoter_nocode_reflection$code, snp_or_cpg_in_promoter_nocode_reflection$interaction)), ]

# perform fisher test
data <- data.frame(c(nrow(snp_or_cpg_in_promoter_code), nrow(snp_or_cpg_in_promoter_nocode)), c(nrow(snp_or_cpg_in_promoter_code_reflection), nrow(snp_or_cpg_in_promoter_nocode_reflection)))
print("Fisher test:")
print(data)
print(unlist(fisher.test(data))[c(4, 1, 2, 3)])
