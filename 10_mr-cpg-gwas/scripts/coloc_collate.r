library(data.table)
library(dplyr)

fn <- paste0("../results/coloc/coloc", 1:191, ".rdata")
l <- list()
for(i in fn)
{
	if(file.exists(i))
	{
		load(i)
		l[[i]] <- res
	} else {
		message(i, " missing")
	}
}

res <- bind_rows(l)
res <- subset(res, !is.na(H4))
ext <- fread("zcat ../results/filtered_gwas_mqtl.txt")
names(ext) <- c("outcome", "snp", "ea", "oa", "eaf", "b", "se", "p", "n")
names(res)[1] <- "snp"

res <- inner_join(res, ext, by=c("snp", "outcome"))

sum(res$p < 1e-10)
sum(res$H4 > 0.8)
sum(res$H4 > 0.8 & res$p < 1e-10)
sig <- res[res$H4 > 0.8 & res$p < 1e-10,]
save(res, file="../results/cpg_trait_coloc.rdata")
