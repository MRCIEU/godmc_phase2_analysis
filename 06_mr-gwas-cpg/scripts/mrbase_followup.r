library(dplyr)
library(TwoSampleMR)
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
			mres <- suppressMessages(mr(dat, metho=c("mr_ivw_radial", "mr_wald_ratio")))
			mres <- subset(mres, ! method %in% "Inverse variance weighted")
			mres$chunk <- chunk
			if(any(! mres$method == "Wald ratio"))
			{
				het <- mr_heterogeneity(dat, method_list=c("mr_ivw_radial")) %>%
					subset(select=c(id.exposure, id.outcome, Q, Q_df, Q_pval))
				mres <- merge(mres, het, by=c("id.exposure", "id.outcome"), all.x=TRUE)
			}
			return(list(mres=mres, dat=dat))
		} else {
			return(NULL)
		}
	} else {
		return(NULL)
	}
}


remove_mhc <- function(x)
{
	require(tidyr)
	x <- tidyr::separate(x, id, c("chr", "pos", "type"), sep=":", remove = FALSE)
	x$pos <- as.numeric(x$pos)
	x <- subset(x, ! (chr == "chr6" & (pos > 25000000 & pos < 35000000)))
	return(x)
}


load("../results/mrbase_sig.rdata")
load("../data/snps_gwas.rdata")


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
	m$all$what <- "all"
	if(nrow(temp2) > 0)
	{
		m$no_mhc <- do_mr(temp2, exp, chunks[i])$mres
		m$no_mhc$what <- "no_mhc"
	}
	if(nrow(temp3) > 0)
	{
		m$mhc <- do_mr(temp3, exp, chunks[i])$mres
		m$mhc$what <- "mhc"
	}
	l[[i]] <- bind_rows(m)
}

res <- bind_rows(l)
res <- subset(res, !is.na(exposure))
res$code <- paste(res$id.exposure, res$id.outcome)
sig$code <- paste(sig$id.exposure, sig$id.outcome)
res <- subset(res, code %in% sig$code)
save(res, file="../results/mrbase_sig_mhc.rdata")



# How many have the same sign and good p-value after excluding mhc
# How many have heterogeneity problems
# How many are not just in MHC


temp <- mutate(res, code=paste(id.exposure, id.outcome)) %>% select(code, b, what)
temp <- spread(temp, key=what, value=b)
key <- subset(res, !duplicated(code), select=c(code, id.exposure, id.outcome, exposure, outcome))
temp <- merge(temp, key, by="code")
cor(temp[,2:4], use="pair")
temp2 <- subset(temp, !is.na(mhc) & !is.na(no_mhc))
sign_agreement <- group_by(temp2, exposure) %>%
	summarise(n=n(), prop_sign=round(sum(sign(mhc) == sign(no_mhc)) / n() * 100)/100) %>% as.data.frame %>% arrange(n)


temp <- mutate(res, code=paste(id.exposure, id.outcome)) %>% select(code, pval, what)
temp <- spread(temp, key=what, value=pval)
key <- subset(res, !duplicated(code), select=c(code, id.exposure, id.outcome, exposure, outcome))
temp <- merge(temp, key, by="code")
cor(temp[,2:4], use="pair")
temp2 <- subset(temp, !is.na(mhc) & !is.na(no_mhc))
sig_agreement <- group_by(temp2, exposure) %>%
	summarise(
		prop_sig_all=sum(all < 1e-5)/n(),
		prop_sig_mhc=sum(mhc < 1e-5)/n(),
		prop_sig_no_mhc=sum(no_mhc < 1e-5)/n()
	) %>% 
	as.data.frame

temp <- mutate(res, code=paste(id.exposure, id.outcome)) %>% select(code, Q_pval, what)
temp <- spread(temp, key=what, value=Q_pval)
key <- subset(res, !duplicated(code), select=c(code, id.exposure, id.outcome, exposure, outcome))
temp <- merge(temp, key, by="code")
cor(temp[,2:4], use="pair")
temp2 <- subset(temp, !is.na(mhc) & !is.na(no_mhc))
qsig_agreement <- group_by(temp2, exposure) %>%
	summarise(n=n(),
		prop_sig_all=sum(all < 0.05)/n(),
		prop_sig_no_mhc=sum(no_mhc < 0.05)/n()
	) %>% 
	as.data.frame %>% arrange(n)



temp <- subset(res, res$what == "no_mhc" & res$Q_pval > 0.05 & res$pval < 1e-7, na.rm=T)
table(temp$exposure) %>% as.data.frame %>% arrange(Freq)