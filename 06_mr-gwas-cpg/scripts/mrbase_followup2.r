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

mr_sign <- function(b_exp, b_out, se_exp=NULL, se_out=NULL, parameters=NULL)
{
	b_exp[b_exp == 0] <- NA
	b_out[b_out == 0] <- NA
	if(sum(!is.na(b_exp) & !is.na(b_out)) < 6)
	{
		return(list(b=NA, se=NA, pval=NA, nsnp=NA))
	}
	x <- sum(sign(b_exp) == sign(b_out), na.rm=TRUE)
	n <- sum(!is.na(b_exp) & !is.na(b_out))

	out <- binom.test(x=x, n=n, p=0.5)
	b <- (out$estimate - 0.5) * 2
	names(b) <- NULL
	pval <- out$p.value
	return(list(b=b, se=NA, pval=pval, nsnp=n))
}
mr_method_list <- function()
{
	a <- list(
		list(
			obj="mr_wald_ratio",
			name="Wald ratio",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		# list(
		# 	obj="mr_meta_fixed_simple",
		# 	name="Fixed effects meta analysis (simple SE)",
		# 	PubmedID="",
		# 	Description="",
		# 	use_by_default=FALSE,
		# 	heterogeneity_test=FALSE
		# ),
		# list(
		# 	obj="mr_meta_fixed",
		# 	name="Fixed effects meta analysis (delta method)",
		# 	PubmedID="",
		# 	Description="",
		# 	use_by_default=FALSE,
		# 	heterogeneity_test=TRUE
		# ),
		# list(
		# 	obj="mr_meta_random",
		# 	name="Random effects meta analysis (delta method)",
		# 	PubmedID="",
		# 	Description="",
		# 	use_by_default=FALSE,
		# 	heterogeneity_test=TRUE
		# ),
		list(
			obj="mr_two_sample_ml",
			name="Maximum likelihood",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_egger_regression",
			name="MR Egger",
			PubmedID="26050253",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_egger_regression_bootstrap",
			name="MR Egger (bootstrap)",
			PubmedID="26050253",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_simple_median",
			name="Simple median",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_weighted_median",
			name="Weighted median",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_penalised_weighted_median",
			name="Penalised weighted median",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_ivw",
			name="Inverse variance weighted",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=TRUE
		),
		list(
			obj="mr_ivw_mre",
			name="Inverse variance weighted (multiplicative random effects)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_ivw_fe",
			name="Inverse variance weighted (fixed effects)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_simple_mode",
			name="Simple mode",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_weighted_mode",
			name="Weighted mode",
			PubmedID="",
			Description="",
			use_by_default=TRUE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_weighted_mode_nome",
			name="Weighted mode (NOME)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_simple_mode_nome",
			name="Simple mode (NOME)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_raps",
			name="Robust adjusted profile score (RAPS)",
			PubmedID="",
			Description="",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		),
		list(
			obj="mr_sign",
			name="Sign concordance test",
			PubmedID="",
			Description="Tests for concordance of signs between exposure and outcome",
			use_by_default=FALSE,
			heterogeneity_test=FALSE
		)
	)
	a <- lapply(a, as.data.frame)
	a <- plyr::rbind.fill(a)
	a <- as.data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
	a$heterogeneity_test <- as.logical(a$heterogeneity_test)
	a$use_by_default <- as.logical(a$use_by_default)
	return(a)
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
			mres <- suppressMessages(mr(dat, metho=c("mr_ivw_radial", "mr_wald_ratio", "mr_sign")))
			mres <- subset(mres, ! method %in% "Inverse variance weighted")
			mres$chunk <- chunk
			if(any(! mres$method == "Wald ratio"))
			{
				het <- mr_heterogeneity(dat, method_list=c("mr_ivw_radial")) %>%
					subset(select=c(id.exposure, id.outcome, Q, Q_df, Q_pval))
				mres <- merge(mres, het, by=c("id.exposure", "id.outcome"), all.x=TRUE)
				outs <- plyr::ddply(dat, .(id.exposure, id.outcome, exposure, outcome), function(x)
				{
					datr <- format_radial_dat(x)
					mod1 <- ivw_radial(datr, alpha=0.05, summary=FALSE)
					nouts <- subset(mod1$data, Qj_Chi > 0.05/nrow(x))$SNP
					datr <- subset(datr, SNP %in% nouts)
					if(nrow(datr) > 1)
					{
						mod2 <- ivw_radial(datr, alpha=0.05, summary=FALSE)
					} else {
						mod2 <- mod1
					}
					data.frame(
						nsnp = c(mod1$df + 1, mod2$df + 1),
						b = c(mod1$coef[1], mod2$coef[1]),
						se = c(mod1$coef[2], mod2$coef[2]),
						pval = c(pnorm(mod1$coef[3], low=FALSE), pnorm(mod2$coef[3], low=FALSE)),
						method = "IVW radial",
						Q = c(mod1$qstatistic, mod2$qstatistic),
						Q_df = c(mod1$df, mod2$df),
						Q_pval = c(pchisq(mod1$qstatistic, mod1$df, low=F), pchisq(mod2$qstatistic, mod2$df, low=F)),
						chunk = chunk,
						what = c("all", "no_outliers")
					)
				})
				mres <- bind_rows(mres, outs)
			}
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

