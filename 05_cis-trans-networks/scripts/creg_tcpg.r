library(dplyr)
library(tidyr)

load("../../results/16/16_clumped.rdata")
# load("../../data/misc/naeem.rdata")
# naeem_keep <- subset(naeem, Flag.discard.keep. == "keep")$probe
# table(unique(clumped$cpg) %in% naeem_keep) / length(unique(clumped$cpg))
# table(naeem_keep %in% naeem$probe) / nrow(naeem)
# clumped <- subset(clumped, cpg %in% naeem_keep)


# Get the trans hits

tcpg <- subset(clumped, !cis & pval < 1e-10)
tcpg <- data_frame(snp=tcpg$snp, tcpg_chr=tcpg$cpgchr, tcpg_pos=tcpg$cpgpos, tcpg_pval=tcpg$pval, tcpg_b=tcpg$Effect, tcpg_chunk=tcpg$chunk, tcpg=tcpg$cpg)

creg <- subset(clumped, cis & pval < 1e-8)
creg <- data_frame(snp=creg$snp, creg_chr=creg$cpgchr, creg_pos=creg$cpgpos, creg_pval=creg$pval, creg_b=creg$Effect, creg_chunk=creg$chunk, creg=creg$cpg)

dat <- inner_join(creg, tcpg, by="snp")
dat <- dplyr::select(dat, creg, tcpg, snp, creg_chr, tcpg_chr, creg_pos, tcpg_pos, creg_pval, tcpg_pval, creg_b, tcpg_b, tcpg_chunk, creg_chunk)
dat$waldratio <- dat$tcpg_b / dat$creg_b
dat$sign <- sign(dat$waldratio)

save(dat, file="../results/creg_tcpg.rdata")
