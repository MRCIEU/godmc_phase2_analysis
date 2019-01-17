library(dplyr)
setwd("/mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/04_conditional_16")


l <- list()
for(i in 1:962)
{
	message(i)
	fn <- paste0("../results_3coh/16_3cohorts/16_", i, "_conditional_cis.rdata")
	if(file.exists(fn))
	{
		load(fn)
		l[[i]] <- clumped
	} else {
		message("missing")
	}
}
conditional.cis <- bind_rows(l)
names(conditional.cis)[names(conditional.cis) == "P-value"] <- "pval"
##
l <- list()
for(i in 1:962)
{
	message(i)
	fn <- paste0("../results_3coh/16_3cohorts/16_", i, "_conditional_trans.rdata")
	if(file.exists(fn))
	{
		load(fn)
		l[[i]] <- clumped
	} else {
		message("missing")
	}
}
conditional.trans <- bind_rows(l)
names(conditional.trans)[names(conditional.trans) == "P-value"] <- "pval"
##

save(conditional.cis, file="../results_3coh/16_3cohorts/16_conditional_3coh_cis.rdata")
save(conditional.trans, file="../results_3coh/16_3cohorts/16_conditional_3coh_trans.rdata")






library(dplyr)
setwd("/mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/04_conditional_16")


l <- list()
for(i in 1:962)
{
	message(i)
	fn <- paste0("../results_3coh/16_3cohorts/16_", i, "_conditional.rdata")
	if(file.exists(fn))
	{
		load(fn)
		l[[i]] <- clumped
	} else {
		message("missing")
	}
}
conditional <- bind_rows(l)
names(conditional)[names(conditional) == "P-value"] <- "pval"

save(conditional, file="../results_3coh/16_3cohorts/16_conditional_3cohorts.rdata")

table(conditional$pJ < 5e-8)
table(conditional$pJ < 5e-14)

cs <- subset(conditional, pJ < 5e-14)

a <- subset(conditional, (cis & pJ < 1e-8) | (!cis & pJ < 1e-14))

b <- group_by(a, cpg, cis) %>%
	dplyr::summarise(n=n())

group_by(b, cis) %>% summarise(mean=mean(n), median=median(n))

group_by(a, cpg) %>% summarise(n=n()) %>% summarise(mean=mean(n), median=median(n))

table(cs$cis)
length(unique(cs$cpg))

range(cs$n)

temp <- conditional[conditional$p < 1e-14 & !conditional$cis,]
length(unique(temp$cpg))
table(table(temp$cpg))

temp <- conditional[conditional$p < 1e-5 & conditional$cis,]
length(unique(temp$cpg))
table(table(temp$cpg))


# Are there any trans with many hits
trans <- subset(conditional, !cis)
x <- names(table(trans$snp)[which.max(table(trans$snp))])
subset(trans, snp==x)$p
bigt <- subset(trans, snp==x)
table(bigt$cpgchr)
bigt %>% as.data.frame
ls()
subset(conditional, snp == x & cis)$p
