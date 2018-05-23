library(dplyr)

ewas <- read.table("../data/EWAS_Catalog_20-02-2018.txt.gz", he=T, sep="\t", stringsAsFactors=FALSE)
ewas$code <- paste(ewas$PMID, ewas$Trait)
b <- group_by(subset(ewas, P < 1e-7), code) %>% summarise(n=n()) %>% arrange(desc(n))

lambda <- function(x)
{
	x <- x[is.finite(x)]
	x[x == 0] <- 1e-300
	ntp <- length(x)
	x <- qchisq(x, 1, lower.tail=FALSE)
	x <- sort(x)
	ppoi <- ppoints(x)
	ppoi <- sort(qchisq(ppoi, df=1, lower.tail=FALSE))
	x <- x[1:ntp]
	ppoi <- ppoi[1:ntp]
	m <- median(x, na.rm=TRUE)/qchisq(0.5, df=1)
	s <- summary( lm(x~0+ppoi) )$coeff
	e <- s[1,1]
	se <- s[1,2]
	return(c(m,e,se))
}

args <- commandArgs(T)
jid <- as.numeric(args[1])
load("../data/snps_gwas.rdata")
load("../results/gwas_summary.rdata")
traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)

grep("Body mass", traits)
jid <- 7
cpglist <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/bmi_cpg.txt", what=character())
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)

min(res$pval)

grep("Cigaret", traits)
jid <- 100
load("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/joehanes.rdata")

load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)


x1 <- data_frame(cpg=joehanes$Probe.ID, eff=joehanes$Effect)
x2 <- data_frame(cpg=res$outcome, b=res$b)
x <- merge(x1, x2)

cor(x$eff, x$b)
summary(lm(eff ~ b, x))



grep("choles", traits)
jid <- 22
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)

cpglist <- subset(ewas, code %in% "28194238 High-density lipoprotein cholesterol")$CpG %>% unique
lambda(subset(res, outcome %in% cpglist)$pval)



jid <- 23
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)

jid <- 24
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)

