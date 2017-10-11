library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

find_bad_snps <- function(x, bim)
{
	notinbim <- subset(x, !snp %in% bim$snp)$snp %>% unique
	x <- subset(x, !snp %in% notinbim)
	index <- match(x$snp, bim$snp)
	bim <- bim[index, ]
	stopifnot(all(bim$snp == x$snp))	

	palin <- (
		(x$Allele1 == "a" & x$Allele2 == "t") |
		(x$Allele1 == "t" & x$Allele2 == "a") |
		(x$Allele1 == "g" & x$Allele2 == "g") |
		(x$Allele1 == "c" & x$Allele2 == "c")) &
		(x$Freq1 > 0.45 & x$Freq1 < 0.55)


	flipped <- (x$Allele1 != bim$a1 & x$Allele2 != bim$a2) &
		(x$Allele2 != bim$a1 & x$Allele1 != bim$a2)

	message(nrow(indel), " indels")
	message(nrow(x), " snps")
	message(sum(palin), " palindromes to exclude")
	message(sum(flipped), " flipped to exclude")

	res <- bim$snp[palin | flipped] %>% unique
	return(res)
}

# load("../data/ref/alleles.rdata")
bim <- fread("../data/ref/eur2.bim")
names(bim) <- c("chr", "snp", "g", "pos", "a1", "a2")
bim$a1 <- tolower(bim$a1)
bim$a2 <- tolower(bim$a2)

bim <- subset(bim, select=c(snp, a1, a2))
save(bim, file="../data/ref/eur2.bim.rdata")

l <- list()
for(i in 1:962)
{
	message(i)
	a <- fread(paste0("zcat ../results/16/16_", i, ".txt.gz"))
	a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
	a$snp2 <- a$snp
	a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
	a$snppos <- as.numeric(a$snppos)
	l[[i]] <- find_bad_snps(a, bim)
}

snplist <- unlist(l) %>% unique
write.table(snplist, file="../data/ref/flipped_snps.rdata", row=F, col=F, qu=F)


