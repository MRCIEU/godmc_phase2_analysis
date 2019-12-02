suppressWarnings(suppressPackageStartupMessages({
	library(gwasvcf)
	library(parallel)
	library(dplyr)
	library(magrittr)
	library(data.table)
	library(ggplot2)
	library(tidyr)
	library(knitr)
	library(pander)
}))
set_bcftools()

load("results/gwas_coloc.rdata")
load("../data/misc/cpg_pos.rdata")
cpgpos$cpgchr <- gsub("chr", "", cpgpos$cpgchr)
res$trait <- tolower(res$trait)
ressig <- subset(res, PP.H4.abf > PP.H3.abf & PP.H4.abf > PP.H2.abf & PP.H4.abf > PP.H1.abf & PP.H4.abf > PP.H0.abf)
ressig <- inner_join(ressig, cpgpos, by="cpg")

traits <- unique(res$trait)

vcfdir <- paste0("~/IGD/data/public/", traits)
file.exists(vcfdir)
vcf <- file.path(vcfdir, paste0(traits, ".vcf.gz"))
stopifnot(all(file.exists(vcf)))

tophits <- mclapply(1:length(vcfdir), function(i)
{
	clump <- scan(file.path(vcfdir[i], "clump.txt"), what="character")
	return(query_gwas(vcf[i], rsid=clump))
}, mc.cores=10)
names(tophits) <- traits

save(tophits, ressig, vcfdir, vcf, traits, cpgpos, file="results/closest.rdata")

