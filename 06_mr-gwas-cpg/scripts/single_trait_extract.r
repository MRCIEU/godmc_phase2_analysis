library(dplyr)

args <- commandArgs(T)
jid <- as.numeric(args[1])
load("../data/snps_gwas.rdata")

traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)
trait <- traits[jid]

l <- list()
for(i in 1:942)
{
	message(i)
	load(paste0("../results/out/out", i, ".rdata"))
	l[[i]] <- subset(res, exposure == trait)
}

res <- bind_rows(l)
save(res, file=paste0("../results/out/gwas", jid, ".rdata"))

