library(dplyr)

load("../results/16/16_clumped.rdata")

assoc <- paste(clumped$cpg, clumped$snp)
temp <- clumped %>% dplyr::select(cpg, snp, cpgchr, cpgpos, cis, snptype)

clumped_aligned <- lapply(1:962, function(i)
{
	message(i)
	load(paste0("../results/16/16_", i, "_aligned.rdata"))
	b <- subset(a, paste(cpg, snp) %in% assoc)
	b$chunk <- i
	names(b)[names(b) == "chr"] <- "snpchr"
	names(b)[names(b) == "pos"] <- "snppos"
	names(b)[names(b) == "Pvalue"] <- "pval"
	b <- inner_join(b, temp, by=c("cpg"="cpg", "snp"="snp"))
	b$PvalueARE <- as.numeric(b$PvalueARE)
	b$PvalueMRE <- as.numeric(b$PvalueMRE)
	return(b)
}) %>% bind_rows() %>% as_tibble()

i <- match(paste(clumped_aligned$snp, clumped_aligned$cpg), paste(clumped$snp, clumped$cpg))
all(clumped_aligned$snp == clumped$snp[i])
table(abs(clumped_aligned$Effect) == abs(clumped$Effect[i]))
table(clumped_aligned$Effect == clumped$Effect[i])


save(clumped_aligned, file="../results/16/16_clumped_aligned.rdata")
