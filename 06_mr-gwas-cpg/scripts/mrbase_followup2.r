library(TwoSampleMR)
library(plyr)
library(dplyr)
library(tidyr)
library(RadialMR)

do_mr <- function(a, b, chunk)
{
	# exposure <- subset(a, exposure == gw)
	exposure <- a
	exposure$SNP <- exposure$id
	# exposure$exposure <- gw
	exposure <- subset(exposure, select=-c(id, trait))

	outcome <- subset(b, snp %in% exposure$SNP)
	if(nrow(outcome) > 0)
	{
		message("Formatting")
		outcome <- format_data(outcome, type="outcome",
			phenotype_col="cpg",
			snp_col="snp",
			beta_col="Effect",
			se_col="StdErr",
			effect_allele_col="Allele1",
			other_allele_col="Allele2",
			eaf_col="Freq1",
			pval_col="P-value",
			samplesize_col="TotalSampleSize"
		)
		outcome$id.outcome <- outcome$outcome
		exposure$SNP <- tolower(exposure$SNP)
		message("Harmonising")
		dat <- suppressMessages(harmonise_data(exposure, outcome, action=2))
		message(nrow(dat))
		if(nrow(dat) > 0)
		{
			message("Analysing")
			res <- suppressMessages(mr(dat, metho=c("mr_ivw", "mr_wald_ratio", "mr_sign", "mr_simple_mode", "mr_weighted_mode", "mr_simple_median", "mr_weighted_median", "mr_egger_regression")))
			het <- suppressMessages(mr_heterogeneity(dat))
			plei <- suppressMessages(mr_pleiotropy_test(dat))
			if(nrow(res) > 0) res$chunk <- chunk
			if(nrow(het) > 0) het$chunk <- chunk
			if(nrow(plei) > 0) plei$chunk <- chunk
			dat$chunk <- chunk
			return(list(res=res, het=het, plei=plei, dat=dat))
		} else {
			return(NULL)
		}
	} else {
		return(NULL)
	}
}


run_chunk <- function(chunk)
{
	res <- list()
	load(paste0("../results/mrbase_dat/dat", chunk, ".rdata"))
	x <- subset(sig, chunk == chunk)
	temp <- subset(a, id.exposure %in% x$id.exposure)

	out <- do_mr(temp, exp, chunk)

	res$res <- out$res
	res$het <- out$het
	res$plei <- out$plei
	dat <- out$dat
	save(dat, file=paste0("../results/mrbase_dat/datf", chunk, ".rdata"))
	return(res)
}


load("../results/mrbase_sig.rdata")
load("../data/snps_gwas.rdata")

out <- list()
chunks <- unique(sig$chunk)
for(i in 1:length(chunks))
{
	message(i)
	out[[i]] <- run_chunk(chunks[i])
}

res <- lapply(out, function(x) x$res) %>% bind_rows %>% as_data_frame
het <- lapply(out, function(x) x$het) %>% bind_rows %>% as_data_frame
plei <- lapply(out, function(x) x$plei) %>% bind_rows %>% as_data_frame

save(res, het, plei, file="../results/mrbase_tophits_full.rdata")


