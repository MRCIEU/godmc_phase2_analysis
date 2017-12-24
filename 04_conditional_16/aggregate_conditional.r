library(dplyr)


l <- list()
for(i in 1:962)
{
	message(i)
	fn <- paste0("../results/16/16_", i, "_conditional.rdata")
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

save(conditional, file="../results/16/16_conditional.rdata")

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
