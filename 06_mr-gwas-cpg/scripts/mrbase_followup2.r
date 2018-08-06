library(plyr)
library(dplyr)
library(TwoSampleMR)
library(tidyr)
library(RadialMR)

remove_mhc <- function(x)
{
	require(tidyr)
	x <- tidyr::separate(x, id, c("chr", "pos", "type"), sep=":", remove = FALSE)
	x$pos <- as.numeric(x$pos)
	x <- subset(x, ! (chr == "chr6" & (pos > 25000000 & pos < 35000000)))
	return(x)
}


format_radial_dat <- function(dat)
{
	dat <- subset(dat, mr_keep)
	format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
}


load("../results/mrbase_sig.rdata")
load("../data/snps_gwas.rdata")

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
			mres <- suppressMessages(mr(dat, metho=c("mr_ivw", "mr_wald_ratio", "mr_sign")))
			# mres <- subset(mres, ! method %in% "Inverse variance weighted")
			mres$chunk <- chunk
			return(list(mres=mres, dat=dat))
		} else {
			return(NULL)
		}
	} else {
		return(NULL)
	}
}



l <- list()
chunks <- unique(sig$chunk)
for(i in 1:length(chunks))
{
	message(i)
	m <- list()
	load(paste0("../results/mrbase_dat/dat", chunks[i], ".rdata"))
	x <- subset(sig, chunk == chunks[i])
	temp <- subset(a, id.exposure %in% x$id.exposure)
	temp2 <- remove_mhc(temp)
	temp3 <- subset(temp, !id %in% temp2$id)
	m$all <- do_mr(temp, exp, chunks[i])$mres
	m$all$what2 <- "all"
	if(nrow(temp2) > 0)
	{
		m$no_mhc <- do_mr(temp2, exp, chunks[i])$mres
		m$no_mhc$what2 <- "no_mhc"
	}
	if(nrow(temp3) > 0)
	{
		m$mhc <- do_mr(temp3, exp, chunks[i])$mres
		m$mhc$what2 <- "mhc"
	}
	l[[i]] <- bind_rows(m)
}

res <- bind_rows(l)
res <- subset(res, !is.na(exposure))
res$code <- paste(res$id.exposure, res$id.outcome)
sig$code <- paste(sig$id.exposure, sig$id.outcome)
res <- subset(res, code %in% sig$code)
save(res, file="../results/mrbase_sig_mhc_sign.rdata")

q()

