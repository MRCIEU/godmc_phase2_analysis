library(dplyr)
library(tidyr)
library(GenomicRanges)

load("../../data/hi-c/nodups.data.Rdata")
load("../../results/16/16_clumped.rdata")
load("../../05_cis-trans-networks/data/snpcontrolsets_selection.rdata")
zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())

clumped <- subset(clumped, cpg %in% zhou & cpgchr != snpchr & pval < 1e-14)

clumped <- subset(clumped, paste(cpg, snp) %in% paste(nodups_data$CpG, nodups_data$SNP))


ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
clumped <- inner_join(clumped, ldinfo, by=c("snp"="SNP"))

save(clumped, file="../data/trans_clumped.rdata")


temp <- subset(clumped, !duplicated(cpg))
cpg <- GRanges(seqnames=temp$cpgchr, ranges=IRanges(temp$cpgpos, temp$cpgpos), strand="*")
names(cpg) <- temp$cpg
temp <- subset(clumped, !duplicated(snp))
snp <- GRanges(seqnames=temp$snpchr, ranges=IRanges(temp$min, temp$max), strand="*")
names(snp) <- temp$snp

save(cpg, snp, file="../data/trans_granges.rdata")

