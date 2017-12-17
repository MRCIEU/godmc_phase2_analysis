library(dplyr)


l <- list()
for(i in 1:962)
{
	message(i)
	fn <- paste0("../results/16/16_", i, "_clumped.rdata")
	if(file.exists(fn))
	{
	load(fn)
	# load(paste0("../results/16/16_", i, "_clumped.rdata"))
	clumped$chunk <- i
	clumped$Pvalue <- as.numeric(clumped$Pvalue)
	clumped$HetPVal <- as.numeric(clumped$HetPVal)
	clumped$PvalueARE <- as.numeric(clumped$PvalueARE)
	clumped$PvalueMRE <- as.numeric(clumped$PvalueMRE)
	l[[i]] <- clumped
	} else {
	message("missing")
	}
}
clumped <- bind_rows(l)
names(clumped)[names(clumped) == "Pvalue"] <- "pval"
save(clumped, file="../results/16/16_clumped.rdata")

# q()

table(clumped$pval < 5e-8)
table(clumped$pval < 5e-14)

sum(clumped$pval < 5e-8 & clumped$cis)
sum(clumped$pval < 5e-14 & !clumped$cis)
q()

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

