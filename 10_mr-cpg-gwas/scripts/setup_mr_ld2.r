# options(dplyr.show_progress = FALSE)

library(dplyr)
library(data.table)

convert_to_chrpos <- function(snp, snp_1kg=snp_1kg)
{
	index <- match(snp, snp_1kg$V2)
	return(snp_1kg$snp[index])
}


message("Loading data")
fn <- read.csv("../../data/gwas/00info.csv")
load("../data/conditional.rdata")
load("../data/snp_1kg.rdata")
extnom <- paste0("../data/extracted/filtered_gwas_mqtl_conditional_", fn$id, ".txt")


for(jid in 1:nrow(fn))
{
	message(jid)
	ext <- fread(extnom[jid])
	names(ext) <- c("snp", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "sample_size.outcome")
	message("subsetting")
	ext <- subset(ext, snp %in% snp_1kg$V2)
	ext$SNP <- convert_to_chrpos(ext$snp, snp_1kg)
	save(ext, file=paste0("../data/extracted/filtered_gwas_mqtl_conditional_ready_", fn$id[jid], ".rdata"))
}


