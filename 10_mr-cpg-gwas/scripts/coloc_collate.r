library(data.table)
library(dplyr)

fn <- paste0("../results/coloc/coloc", 1:164, ".rdata")
l <- list()
for(i in fn)
{
	load(i)
	l[[i]] <- res
}

res <- bind_rows(l)
res <- subset(res, !is.na(H4))
ext <- fread("zcat ../results/filtered_gwas_mqtl.txt")
names(ext) <- c("outcome", "snp", "ea", "oa", "eaf", "b", "se", "p", "n")
names(res)[1] <- "snp"

res <- inner_join(res, ext, by=c("snp", "outcome"))
save(res, file="../results/cpg-trait.rdata")
