library(tidyverse)
load("../data/annotations.rdata")
load("../data/trans_clumped.rdata")


annokeep <- intersect(unique(cpgres$anno), unique(snpres$anno))
cpgres <- subset(cpgres, anno %in% annokeep)
snpres <- subset(snpres, anno %in% annokeep)
snpres <- subset(snpres, !duplicated(paste(anno, snp)))
cpgres <- subset(cpgres, !duplicated(paste(anno, cpg)))

cpgcount <- group_by(cpgres, anno) %>%
	summarise(ncpg=n())

snpcount <- group_by(snpres, anno) %>%
	summarise(nsnp=n())

counts <- inner_join(cpgcount, snpcount)

counts$pcpg <- counts$ncpg / nrow(clumped)
counts$psnp <- counts$nsnp / nrow(clumped)

expected <- counts$psnp %*% t(counts$pcpg) 

expected[1:10, 1:10] * nrow(clumped)
mat[1:10,1:10]

expected <- expected %>%
	reshape2::melt()

load("../results/matrix/m0.rdata")
res <- reshape2::melt(mat)
res$expected <- expected$value

res$pval1 <- pbinom(res$value, nrow(clumped), res$expected)
res$pval2 <- pbinom(res$value, nrow(clumped), res$expected, low=FALSE)
res$pval <- pmin(res$pval1, res$pval2)


load("../results/matrix/m1.rdata")
res_perm <- reshape2::melt(mat)
res_perm$expected <- expected$value

res_perm$pval1 <- pbinom(res_perm$value, nrow(clumped), res_perm$expected)
res_perm$pval2 <- pbinom(res_perm$value, nrow(clumped), res_perm$expected, low=FALSE)
res_perm$pval <- pmin(res_perm$pval1, res_perm$pval2)

sum(res_perm$pval < (0.05/nrow(res_perm)))
sum(res$pval < (0.05/nrow(res)))


load("../results/matrix/m2.rdata")
res_perm <- reshape2::melt(mat)
res_perm$expected <- expected$value

res_perm$pval1 <- pbinom(res_perm$value, nrow(clumped), res_perm$expected)
res_perm$pval2 <- pbinom(res_perm$value, nrow(clumped), res_perm$expected, low=FALSE)
res_perm$pval <- pmin(res_perm$pval1, res_perm$pval2)

sum(res_perm$pval < (0.05/nrow(res_perm)))
sum(res$pval < (0.05/nrow(res)))


