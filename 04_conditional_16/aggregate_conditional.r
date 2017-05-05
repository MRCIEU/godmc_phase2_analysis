library(dplyr)


l <- list()
for(i in 1:962)
{
	load(paste0("../results/16/16_", i, "_conditional.rdata"))
	l[[i]] <- clumped
}
conditional <- bind_rows(l)
names(conditional)[names(conditional) == "P-value"] <- "pval"

save(conditional, file="../results/16/16_conditional.rdata")

table(conditional$pval < 5e-8)
table(conditional$pval < 5e-14)

cs <- subset(conditional, pval < 5e-14)

table(cs$cis)
length(unique(cs$cpg))

range(cs$TotalSampleSize)

temp <- conditional[conditional$pval < 1e-14 & !conditional$cis,]
length(unique(temp$cpg))
table(table(temp$cpg))

temp <- conditional[conditional$pval < 1e-5 & conditional$cis,]
length(unique(temp$cpg))
table(table(temp$cpg))


# Are there any trans with many hits
trans <- subset(conditional, !cis)
x <- names(table(trans$snp)[which.max(table(trans$snp))])
subset(trans, snp==x)$pval
bigt <- subset(trans, snp==x)
table(bigt$cpgchr)
bigt %>% as.data.frame
ls()
subset(conditional, snp == x & cis)$pval
