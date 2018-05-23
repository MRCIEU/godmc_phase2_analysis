# How many have the same sign and good p-value after excluding mhc
# How many have heterogeneity problems
# How many are not just in MHC


# Criteria
# No heterogeneity in IVW
# 


library(plyr)
library(dplyr)
library(TwoSampleMR)
library(tidyr)
library(RadialMR)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- subset(ann, select=c(chr, pos))
ann$outcome <- rownames(ann)

# load("../results/mrbase_sig.rdata")
load("../data/snps_gwas.rdata")
load("../results/mrbase_sig_mhc.rdata")
load("../results/mrbase_sig.rdata")

sig$code <- paste(sig$id.exposure, sig$id.outcome)

a0 <- filter(a, data_source.exposure == "mrbase") %>% group_by(id.exposure, exposure) %>% summarise(n=n()) %>% arrange(exposure, n)

# Some traits are duplicated, look by eye to remove the out of date ones
data.frame(substr(as.character(a1$exposure), start=1, stop=50), a1$n)
remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)
a1 <- a0[-remove_ids,]

keep_traits <- a1$exposure
res <- subset(res, exposure %in% keep_traits)
sig <- subset(sig, exposure %in% keep_traits)
res2 <- subset(res, pval < 1.4e-7)
t0 <- subset(res2, !duplicated(exposure))

sigcodes <- subset(sig, pval < 1.4e-7)$code %>% unique
sigcodes_single <- subset(sig, nsnp == 1 & pval < 1.4e-7)$code %>% unique
sigcodes_multiple <- subset(sig, nsnp > 1 & pval < 1.4e-7)$code %>% unique



# The single instrument cases
# Remove MHC examples


t1 <- subset(res, code %in% sigcodes_single & is.na(what))
table(t1$exposure, t1$what2) %>% as.data.frame %>% spread(Var2, Freq) %>% arrange(all)

subset(t1, grepl("Intracranial volume", exposure))$outcome %>% unique

# Let's extract the original data for these associations
a_single <- subset(a, id.exposure %in% t1$id.exposure)

dat <- group_by(t1, chunk) %>%
do({
	x <- .
	load(paste0("../results/mrbase_dat/dat", x$chunk[1], ".rdata"))
	exp <- subset(exp, snp %in% a_single$id)
	outcome <- suppressMessages(format_data(exp, type="outcome",
		phenotype_col="cpg",
		snp_col="snp",
		beta_col="Effect",
		se_col="StdErr",
		effect_allele_col="Allele1",
		other_allele_col="Allele2",
		eaf_col="Freq1",
		pval_col="P-value",
		samplesize_col="TotalSampleSize"
	))
	outcome$id.outcome <- outcome$outcome
	exposure <- a_single
	exposure$SNP <- tolower(exposure$id)
	outcome$SNP <- tolower(outcome$SNP)
	dat <- suppressMessages(harmonise_data(exposure, outcome))
	dat$code <- paste(dat$id.exposure, dat$id.outcome)
	out <- subset(dat, code %in% sigcodes_single)
	out
})

dat <- subset(dat, !duplicated(dat$code))
datcont <- subset(dat, units.exposure != "log odds")
datbin <- subset(dat, units.exposure == "log odds")
datbin$rsq.exposure <- get_r_from_lor(datbin$beta.exposure, datbin$eaf.outcome, datbin$ncase.exposure, datbin$ncontrol.exposure, 0.01)
datbin$rsq.outcome <- 2 * datbin$beta.outcome^2 * datbin$eaf.outcome * (1 - datbin$eaf.outcome)
# datcont$rsq.exposure <- 2 * datcont$beta.exposure^2 * datcont$eaf.outcome * (1 - datcont$eaf.outcome)
datcont$rsq.exposure <- get_r_from_pn(datcont$pval.exposure, datcont$samplesize.exposure)^2
datcont$rsq.outcome <- 2 * datcont$beta.outcome^2 * datcont$eaf.outcome * (1 - datcont$eaf.outcome)

dat <- rbind(datbin, datcont)
group_by(dat, exposure) %>%
summarise(n=n(), nforward=sum(correct_dir & steiger_p < 0.05), nrev=sum(!correct_dir & steiger_p < 0.05)) %>% as.data.frame

dat$correct_dir <- dat$rsq.exposure > dat$rsq.outcome
dat$steiger_p <- psych::r.test(n = dat$samplesize.exposure, n2 = dat$samplesize.outcome, r12 = dat$rsq.exposure, r34 = dat$rsq.outcome)$p

dats <- subset(dat, select=c(code, SNP, effect_allele.exposure, other_allele.exposure, correct_dir, steiger_p))
a_temp <- subset(a, select=c(SNP, id))
a_temp$id <- tolower(a_temp$id)
a_temp <- tidyr::separate(a_temp, id, c("chr", "pos", "type"), sep=":", remove = FALSE)
a_temp$pos <- as.numeric(a_temp$pos)
a_temp$MHC <- a_temp$chr == "chr6" & (a_temp$pos > 25000000 & a_temp$pos < 35000000)
single_tab <- merge(subset(sig, select=c(code, exposure, outcome, b, se, pval, nsnp)), dats, by="code") 
single_tab <- merge(a_temp, single_tab, by.x="id", by.y="SNP") %>% subset(select=c(exposure, outcome, SNP, chr, pos, MHC, effect_allele.exposure, other_allele.exposure, b, se, pval, correct_dir, steiger_p))


names(single_tab) <- c("exposure", "outcome", "SNP", "chr", "pos", "mhc", "A1", "A2", "Wald ratio", "SE", "p", "correct_dir", "steiger_p")

single_tab <- single_tab[order(single_tab$mhc, single_tab$exposure, single_tab$p), ]

save(single_tab, file="../results/mrbase_single_instrument_results.rdata")



## Directionality



# Significant associations with no outliers
# Any outliers?
# 




sigcodes_multiple2 <- subset(res, what=="all" & what2 == "all" & pval < 1.4e-7)$code %>% unique

t2 <- subset(res, code %in% sigcodes_multiple2) %>% arrange(code)
table(t2$exposure, t2$what2) %>% as.data.frame %>% spread(Var2, Freq) %>% arrange(all)
head(t2, 10)
t2$what <- as.character(t2$what)
t2$what[is.na(t2$what)] <- "na"

x <- subset(t2, code == "3cI0qp cg17713376")

t2 <- group_by(t2, code) %>%
do({
	x <- .

	x <- subset(x, !is.na(se) & what != "no_outliers")
	x <- subset(x, !(what == "na" & method == "IVW radial"))
	x$what <- "all"
	x
})

# Start with MHC - are there any that are completely lost when MHC is removed

# Which traits have only MHC SNPs?
t4 <- group_by(t2, code) %>%
summarise(n=n(), what=paste(what2, collapse=" "))
t4_mhc <- subset(t2, code %in% subset(t4, what=="all mhc")$code) %>% arrange(code)
mhc_code <- t4_mhc$code %>% unique

# Which traits have MHC and non-MHC?
t4_both <- subset(t2, code %in% subset(t4, what=="all no_mhc mhc")$code) %>% arrange(code)

# What associations are significantly different after removing MHC?
dwh <- function(b,se)
{
	h <- (b[1]-b[2])^2 / (se[1]^2-se[2]^2)
	return(pchisq(h, 1, low=F))
}

t4_dwh <- filter(t4_both, what2 != "all") %>%
	group_by(code) %>%
	summarise(
		n = n(), 
		dwh_p = dwh(b,se), 
		diff = dwh_p < 0.05,
		sign_agreement = sign(b[1]) == sign(b[2]),
		pval_thresh = sum(pval < 0.05),
		fi = what2[1],
		pval_no_mhc = pval[1] < 0.05,
		pval_mhc = pval[2] < 0.05
	)
t4_dwh$fdr <- p.adjust(t4_dwh$dwh_p, "fdr") < 0.05


with(t4_dwh, table(paste(pval_no_mhc, pval_mhc), sign_agreement))


t4_both <- merge(t4_both, t4_dwh, by="code")
group_by(t4_both, exposure) %>%
	summarise(n=n(), nfail=sum(diff), signag = sum(sign_agreement)/n()) %>% as.data.frame %>% arrange(signag)


# Traits with no MHC instruments

t5 <- subset(res, code %in% sigcodes_multiple2 & !code %in% mhc_code & !code %in% t4_dwh$code & !is.na(what) & what2 == "all") %>% arrange(code)


## Make table
# Traits with no MHC
# Traits with both MHC and non-MHC
# Traits with only MHC

t5$section <- "No MHC instruments"
t4_both$section <- "MHC and non-MHC instruments"
t4_mhc$section <- "MHC instruments only"
multiple_tab <- bind_rows(subset(t5, what == "all"), t4_both, subset(t4_mhc, what2=="mhc"))
multiple_tab <- subset(multiple_tab, select=c(section, exposure, outcome, nsnp, b, se, pval, Q, Q_pval,what2))
save(multiple_tab, file="../results/mrbase_multiple_instrument_results.rdata")










library(missMethyl)
library(GO.db)
library(KEGGREST)

l <- subset(t5, grepl("Small ", exposure))$outcome %>% unique
go.res2 <- gometh(sig.cpg=l, array.type=c("450K"),collection="KEGG")
go.res2 <- gometh(sig.cpg=cpgs[,1], array.type=c("450K"),collection="KEGG")
go.res2[go.res2$FDR<0.9,]





t2s <- subset(t2, !is.na(what) & what2 == "all")
t2_all <- subset(t2s, what == "all")
table(t2_all$Q_pval < 0.05)


t3 <- subset(res, pval < 1.1e-9)
t4 <- subset(t3, !duplicated(paste(exposure, method)))

s1 <- subset(sig, pval < 1.1e-9)
s2 <- subset(s1, method == "Wald ratio" & !duplicated(exposure))
s3 <- subset(s1, method != "Wald ratio" & !duplicated(exposure))



cigs <- subset(t1, grepl("Cigarett", exposure) & what2=="all")



temp <- filter(res, method == "IVW radial", what2=="all", !is.na(what)) %>%
	group_by(code) %>%
	summarise(outs = Q_df[1] != Q_df[2])

temp2 <- subset(res, !is.na(what) & what2 == "all" & code %in% subset(temp, outs)$code) %>%
	group_by(code) %>%
	summarise(pval1 = pval[1] < 1e-7, pval2 = pval[2] < 1e-7, pval_less = pval2 < pval1, sign = sign(b[1]) == sign(b[2]), conc = b[2] > b[1] - se[1] * 1.96 & b[2] < b[1] + se[1] * 1.96)






sig$code <- paste(sig$id.exposure, sig$id.outcome)
thresh <- 0.05/(350000*140)
sig2 <- subset(sig, pval < thresh)
res <- subset(res, code %in% sig2$code)

temp <- mutate(res, code=paste(id.exposure, id.outcome)) %>% select(code, b, what)
temp <- spread(temp, key=what, value=b)
key <- subset(res, !duplicated(code), select=c(code, id.exposure, id.outcome, exposure, outcome))
temp <- merge(temp, key, by="code")
cor(temp[,2:4], use="pair")
temp2 <- subset(temp, !is.na(mhc) & !is.na(no_mhc))
sign_agreement <- group_by(temp2, exposure) %>%
	summarise(n=n(), prop_sign=round(sum(sign(mhc) == sign(no_mhc)) / n() * 100)/100) %>% as.data.frame %>% arrange(n)
sum(sign_agreement$n * sign_agreement$prop_sign) / sum(sign_agreement$n)

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
