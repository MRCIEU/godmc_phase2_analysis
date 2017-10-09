library(dplyr)

# setwd("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/scripts/godmc_phase2_analysis/06_mr-gwas-cpg/results")
nom <- paste0("../results/out/", list.files("../results/out", pattern="*.rdata$"))

l <- list()
for(i in 1:length(nom))
{
	load(nom[i])
	message(i)
	message("nrow: ", nrow(res))
	message("ncpg: ", length(unique(res$outcome)))
	message("nexp: ", length(unique(res$exposure)))
	message("nsig: ", sum(res$pval < 1e-5))
	if(nrow(res) > 0)
	{
		res$nom <- nom[i]
		l[[i]] <- subset(res, pval < 1e-5, select=c(exposure, outcome, b, se, pval, nsnp))
	}
}

res <- bind_rows(l)
save(res,file="~/repo/godmc_phase2_analysis/06_mr-gwas-cpg/results/tophits.rdata")

length(unique(res$exposure))
length(unique(res$outcome))

sum(res$pval < 1e-10)

a <- subset(res, pval < 1e-7)

table(table(a$exposure))
unique(a$exposure)

length(unique(a$outcome))

subset(a, grepl("Urate", exposure))


