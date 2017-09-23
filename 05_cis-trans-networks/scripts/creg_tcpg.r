library(dplyr)
library(tidyr)
library(doParallel)

load("../../results/16/16_clumped.rdata")
# load("../../data/misc/naeem.rdata")
# naeem_keep <- subset(naeem, Flag.discard.keep. == "keep")$probe
# table(unique(clumped$cpg) %in% naeem_keep) / length(unique(clumped$cpg))
# table(naeem_keep %in% naeem$probe) / nrow(naeem)
# clumped <- subset(clumped, cpg %in% naeem_keep)


# Get the trans hits

tcpg <- subset(clumped, !cis & pval < 1e-10)
tcpg <- data_frame(snp=tcpg$snp, tcpg_chr=tcpg$cpgchr, tcpg_pos=tcpg$cpgpos, tcpg_pval=tcpg$pval, tcpg_b=tcpg$Effect, tcpg_chunk=tcpg$chunk, tcpg=tcpg$cpg)

# creg <- subset(clumped, cis & pval < 1e-8)
# creg <- data_frame(snp=creg$snp, creg_chr=creg$cpgchr, creg_pos=creg$cpgpos, creg_pval=creg$pval, creg_b=creg$Effect, creg_chunk=creg$chunk, creg=creg$cpg)

cis_radius <- 1000000
(no_cores <- detectCores() - 1)
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores, type="FORK")
result <- parLapply(cl, 1:962, function(i)
{
	message(i)
	load("../../03_clumping_16/cpg_pos.rdata")
	a <- fread(paste0("zcat ../../results/16/16_", i, ".txt.gz"))
	a$Pvalue <- as.numeric(a$Pvalue)
	a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
	a$snp2 <- a$snp
	a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
	a$snppos <- as.numeric(a$snppos)
	a <- inner_join(a, cpgpos, by=c("cpg"))
	a$cis <- FALSE
	a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE
	b <- subset(a, cis & Pvalue < 0.05/ (nrow(tcpg) * 10) & snp %in% tcpg$snp)
	return(data_frame(snp=b$snp, creg_chr=b$cpgchr, creg_pos=b$cpgpos, creg_pval=b$Pvalue, creg_b=b$Effect, creg_chunk=i, creg=b$cpg))

})
stopCluster(cl)

creg <- bind_rows(result)

table(unique(tcpg$snp) %in% creg$snp)

dat <- inner_join(creg, tcpg, by="snp")
dat <- dplyr::select(dat, creg, tcpg, snp, creg_chr, tcpg_chr, creg_pos, tcpg_pos, creg_pval, tcpg_pval, creg_b, tcpg_b, tcpg_chunk, creg_chunk)
dat$waldratio <- dat$tcpg_b / dat$creg_b
dat$sign <- sign(dat$waldratio)

save(dat, file="../results/creg_tcpg.rdata")
