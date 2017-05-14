library(dplyr)


l <- list()
for(i in 1:962)
{
	message(i)
	load(paste0("../results/16/16_", i, "_clumped.rdata"))
	l[[i]] <- clumped
}
clumped <- bind_rows(l)
names(clumped)[names(clumped) == "Pvalue"] <- "pval"

save(clumped, file="../results/16/16_clumped.rdata")

table(clumped$pval < 5e-8)
table(clumped$pval < 5e-14)

cs <- subset(clumped, pval < 5e-14)

table(cs$cis)
length(unique(cs$cpg))

range(cs$TotalSampleSize)

temp <- clumped[clumped$pval < 1e-5 & clumped$cis,]
length(unique(temp$cpg))

table(table(clumped$cpg))


# Are there any trans with many hits
trans <- subset(clumped, !cis)
x <- names(table(trans$snp)[which.max(table(trans$snp))])
subset(trans, snp==x)$pval
bigt <- subset(trans, snp==x)
table(bigt$cpgchr)
bigt %>% as.data.frame
ls()
subset(clumped, snp == x & cis)$pval
