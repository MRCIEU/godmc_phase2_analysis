library(dplyr)

args <- commandArgs(T)
jid <- as.numeric(args[1])
load("../data/snps_gwas.rdata")

traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)
# trait <- traits[jid]

l <- list()
for(i in 1:length(traits))
{
	l[[i]] <- list()
}
for(i in 1:942)
{
	message(i)
	load(paste0("../results/out/out", i, ".rdata"))
	for(j in 1:length(traits))
	{
		cat(".")
		l[[j]][[i]] <- subset(res, exposure == traits[j])
	}
	cat("\n")
}

for(i in 1:length(traits))
{
	message(i)
	res <- bind_rows(l[[i]])
	save(res, file=paste0("../results/out/gwas", i, ".rdata"))
}

