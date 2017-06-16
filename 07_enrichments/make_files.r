library(tidyverse)

load("../results/enrichments/clumpedwithldregion.rdata")
load("../results/enrichments/cis_trans.rdata")
load("../results/enrichments/mqtl_pchic.rdata")
load("../results/16/16_clumped.rdata")


dat <- data_frame(
	cpgchr=cl.out$cpgchr,
	cpgpos=cl.out$cpgpos,
	cpgname=cl.out$cpg,
	snpchr=cl.out$snpchr,
	snppos=cl.out$snppos,
	snpname=cl.out$snp,
	effect_allele=cl.out$Allele1,
	other_allele=cl.out$Allele2,
	effect_allele_freq=cl.out$Freq1,
	mqtl_effect=cl.out$Effect,
	mqtl_pval=cl.out$pval,
	meta_directions=cl.out$Direction,
	Isq=cl.out$HetISq,
	samplesize=cl.out$TotalSampleSize,
	cis=cl.out$cis,
	ld80_start=cl.out$start_bp,
	ld80_end=cl.out$stop_bp,
	ld80_proxies=cl.out$nproxies,
	ld80_dist=ld80_end-ld80_start
)

dat$ld80_proxies[is.na(dat$ld80_proxies)] <- 0

cpg_clusters <- inner_join(
	select(dat, cpgchr, cpgpos, cpgname) %>% filter(!duplicated(cpgname)),
	m,
	by=c("cpgname"="cpg")
)


dat$code <- paste(dat$cpgname, dat$snpname)
names(abc) <- c("code", "pchic")
dat <- left_join(dat, abc, by="code")

save(dat, cpg_clusters, file="../results/enrichments/for_christoph.rdata")