library(dplyr)
library(igraph)
library(coloc)
library(data.table)
library(tidyr)
library(TwoSampleMR)

run_coloc <- function(dat, type1, type2)
{

	# H0 (no causal variant), 
	# H1 (causal variant for trait 1 only), 
	# H2 (causal variant for trait 2 only), 
	# H3 (two distinct causal variants)
	# H4 (one common causal variant)

	require(coloc)
	dataset1 <- list(
		pvalues=dat$pval.exposure,
		N=dat$samplesize.exposure,
		beta=dat$beta.exposure,
		varbeta=dat$se.exposure^2,
		snp=dat$SNP,
		type=type1,
		MAF=dat$eaf.exposure
	)
	dataset2 <- list(
		pvalues=dat$pval.outcome,
		N=dat$samplesize.outcome,
		beta=dat$beta.outcome,
		varbeta=dat$se.outcome^2,
		snp=dat$SNP,
		type=type2,
		MAF=dat$eaf.outcome
	)
	dataset2$MAF[is.na(dataset2$MAF)] <- dataset1$MAF[is.na(dataset2$MAF)]
	a <- suppressMessages(coloc.abf(dataset1, dataset2))
	return(data.frame(
		exposure=dat$exposure[1],
		outcome=dat$outcome[1],
		nsnp=a$summary[1],
		H0=a$summary[2],
		H1=a$summary[3],
		H2=a$summary[4],
		H3=a$summary[5],
		H4=a$summary[6]
	))
}

format_metal <- function(x, type)
{
	x <- suppressMessages(format_data(
		x,
		type=type,
		snp_col="snp",
		phenotype_col="cpg",
		effect_allele_col="Allele1",
		other_allele_col="Allele2",
		eaf_col="Freq1",
		beta_col="Effect",
		se_col="StdErr",
		pval_col="Pvalue",
		samplesize_col="TotalSampleSize"
	))
	return(x)
}

coloc_wrapper <- function(creg, tcpg, creg_chunk, tcpg_chunk)
{
	creg_chunk <- paste0("zcat ../../results/16/16_", creg_chunk, ".txt.gz")
	tcpg_chunk <- paste0("zcat ../../results/16/16_", tcpg_chunk, ".txt.gz")
	creg_dat <- fread(creg_chunk) %>%
		separate(MarkerName, c("snp", "cpg"), sep="_") %>%
		filter(cpg == creg) %>%
		format_metal(type="exposure")

	tcpg_dat <- fread(tcpg_chunk) %>%
		separate(MarkerName, c("snp", "cpg"), sep="_") %>%
		filter(cpg == tcpg, tolower(snp) %in% creg_dat$SNP) %>%
		format_metal(type="outcome")

	dat <- try(harmonise_data(creg_dat, tcpg_dat, action=1))
	if(class(dat) == 'try-error')
	{
		message("NOT RUNNING")
		res <- NULL
	} else {
		res <- run_coloc(dat, "quant", "quant")
	}
	return(res)
}


load("../../results/16/16_clumped.rdata")
load("../results/creg_tcpg.rdata")

args <- commandArgs(T)
jid <- as.numeric(args[1])
num <- as.numeric(args[2])

out <- paste0("../results/coloc", jid, ".rdata")
if(file.exists(out)) q()

first <- (jid - 1) * num + 1
last <- min(jid * num, nrow(dat))
dat <- dat[first:last,]

res <- list()
for(i in 1:nrow(dat))
{
	message(i)
	res[[i]] <- with(dat[i,], coloc_wrapper(creg, tcpg, creg_chunk, tcpg_chunk))
}

res <- bind_rows(res)
save(res, file=out)
