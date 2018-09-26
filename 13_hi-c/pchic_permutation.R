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

# baitName contains NA's, so just remove
pchic <- pchic[, !colnames(pchic) == "baitName"] 

# make unique snp cpg code
clumped$code <- paste(clumped$snp, clumped$cpg)

# make granges for cpgs locations and snp locations
clumped.cpg <- GRanges(seqnames = clumped$cpgchr, 
                       ranges = IRanges(start = clumped$cpgpos, end = clumped$cpgpos, names = clumped$cpg))
clumped.snp <- GRanges(seqnames = clumped$snpchr, 
                       ranges = IRanges(start = clumped$snppos, end = clumped$snppos, names = clumped$snp))

# make unique pchic code
pchic$interaction <- paste(pchic$baitID, pchic$oeStart)

# calculate interaction region distances relative to start of promoter
pchic$oeStartDist <- pchic$baitStart - pchic$oeStart
pchic$oeEndDist <- pchic$baitStart - pchic$oeEnd

# calculate overlap
overlap <- function(pchic) {
  
  # make granges for promoters and interaction regions
  pchic.bait <- GRanges(seqnames = paste0("chr", pchic$baitChr), 
                        ranges = IRanges(start = pchic$baitStart, end = pchic$baitEnd))
  values(pchic.bait) <- pchic
  pchic.oe <- GRanges(seqnames = paste0("chr", pchic$oeChr), 
                      ranges = IRanges(start = pchic$oeStart, end = pchic$oeEnd))
  values(pchic.oe) <- pchic
  
  # calculate overlap between promoter and cpg
  overlaps <- findOverlaps(pchic.bait, clumped.cpg, ignore.strand = T)
  bait.cpg <- values(pchic.bait[queryHits(overlaps)])
  bait.cpg$CpG <- names(clumped.cpg[subjectHits(overlaps)])
  
  # calculate overlap between promoter and snp
  overlaps <- findOverlaps(pchic.bait, clumped.snp, ignore.strand = T)
  bait.snp <- values(pchic.bait[queryHits(overlaps)])
  bait.snp$SNP <- names(clumped.snp[subjectHits(overlaps)])
  
  # calculate overlap between interacting region and cpg
  overlaps <- findOverlaps(pchic.oe, clumped.cpg, ignore.strand = T)
  oe.cpg <- values(pchic.oe[queryHits(overlaps)])
  oe.cpg$CpG <- names(clumped.cpg[subjectHits(overlaps)])
  
  # calculate overlap between interacting region and snp
  overlaps <- findOverlaps(pchic.oe, clumped.snp, ignore.strand = T)
  oe.snp <- values(pchic.oe[queryHits(overlaps)])
  oe.snp$SNP <- names(clumped.snp[subjectHits(overlaps)])
  
  # merge snp in promoter and cpg in interacting region
  snp.in.promoter <- oe.cpg[match(bait.snp$interaction, oe.cpg$interaction), ]
  snp.in.promoter$SNP <- bait.snp$SNP
  snp.in.promoter$code <- paste(snp.in.promoter$SNP, snp.in.promoter$CpG)

  # merge cpg in promoter and snp in interacting region
  cpg.in.promoter <- oe.snp[match(bait.cpg$interaction, oe.snp$interaction), ]
  cpg.in.promoter$CpG <- bait.cpg$CpG
  cpg.in.promoter$code <- paste(cpg.in.promoter$SNP, cpg.in.promoter$CpG)
  
  # combine and make sure snp-cpg pair is in clumped and no duplicates in combination of snp-cpg and promoter-interacting_region
  res <- rbind(snp.in.promoter, cpg.in.promoter[, c(1:32, 34, 33, 35)])
  res <- res[res$code %in% clumped$code, ]
  na.omit(res[!duplicated(paste0(res$interaction, res$code)), ])
  
}

# only consider interactions where either snp in promoter or cpg in promoter
pchic.bait <- GRanges(seqnames = paste0("chr", pchic$baitChr), 
                      ranges = IRanges(start = pchic$baitStart, end = pchic$baitEnd))
values(pchic.bait) <- pchic

overlaps <- findOverlaps(pchic.bait, clumped.cpg, ignore.strand = T)
bait.cpg <- values(pchic.bait[queryHits(overlaps)])

overlaps <- findOverlaps(pchic.bait, clumped.snp, ignore.strand = T)
bait.snp <- values(pchic.bait[queryHits(overlaps)])

pchic.promoter <- rbind(bait.cpg, bait.snp)
rm(pchic.bait, overlaps, bait.cpg, bait.snp)  
gc()

pchic.promoter <- pchic.promoter[!duplicated(pchic.promoter$interaction), ]

enrich <- bplapply(1:1000, function(i) {
  
  # swap IDs
  swaps <- lapply(unique(pchic.promoter$baitChr), function(chr) {
    
    pchic.promoter <- pchic.promoter[pchic.promoter$baitChr == chr, ]
    pchic.promoter <- pchic.promoter[!duplicated(pchic.promoter$baitID), ]
    data.frame(baitID = pchic.promoter$baitID, baitStart=pchic.promoter$baitStart, baitEnd=pchic.promoter$baitEnd, swap = sample(pchic.promoter$baitID, length(pchic.promoter$baitID)))
    
  })
  swaps <- do.call(rbind, swaps)

  # for each promoter get new interacting regions from swapped ID
  pchic.swapped <- lapply(swaps$baitID, function(baitID) {
    
    swap <- swaps$swap[na.omit(swaps$baitID == baitID)]
    pchic.promoter <- pchic.promoter[na.omit(pchic.promoter$baitID == swap), ]
    pchic.promoter$baitID <- baitID
    pchic.promoter$baitStart <- swaps$baitStart[na.omit(swaps$baitID == baitID)]
    pchic.promoter$baitEnd <- swaps$baitEnd[na.omit(swaps$baitID == baitID)]
    pchic.promoter
    
  })
  pchic.swapped <- do.call(rbind, pchic.swapped)
  pchic.swapped$oeStart <- pchic.swapped$baitStart - pchic.swapped$oeStartDist
  pchic.swapped$oeEnd <- pchic.swapped$baitStart - pchic.swapped$oeEndDist
  pchic.swapped$interaction <- paste(pchic.swapped$baitID, pchic.swapped$oeStart)
  
  #remove interactions in swapped data where interacting region outside chromosome
  pchic.swapped <- lapply(unique(pchic.swapped$baitChr), function(chr) {
    
    pchic.swapped <- pchic.swapped[pchic.swapped$baitChr == chr, ]
    pchic.swapped$keep <- 1
    pchic.swapped$keep[pchic.swapped$oeStart < min(pchic$oeStart) | pchic.swapped$oeStart > max(pchic$oeStart)] <- NA
    pchic.swapped$keep[pchic.swapped$oeEnd < min(pchic$oeEnd) | pchic.swapped$oeEnd > max(pchic$oeEnd)] <- NA
    pchic.swapped[!is.na(pchic.swapped$keep), ]
    
  })
  pchic.swapped <- do.call(rbind, pchic.swapped)

  #get overlap for both original and swapped data
  overlap.pchic <- overlap(pchic.promoter[, 1:32])
  overlap.pchic.swapped <- overlap(pchic.swapped[, 1:32])
  overlap.pchic <- data.frame(total = nrow(pchic.promoter), enrich = nrow(overlap.pchic))
  overlap.pchic.swapped <- data.frame(total = nrow(pchic.swapped), enrich = nrow(overlap.pchic.swapped))
  data.frame(overlap.pchic, overlap.pchic.swapped)
  
}, BPPARAM=MulticoreParam(50))

#calculate permutation p-value
enrich <- do.call(rbind, enrich)
enrich2 <- enrich[, 2] / enrich[, 1] <= enrich[, 4] / enrich[, 3]
print(paste0("p = ", sum(enrich2) / length(enrich2)))
