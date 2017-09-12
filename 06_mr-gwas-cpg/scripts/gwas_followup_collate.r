library(dplyr)

setwd("../results")
nom <- list.files(pattern="*.rdata$")
nom <- nom[grepl("fol", nom)]

l <- list()
for(i in 1:length(nom))
{
	load(nom[i])
	message(i)
	message("nrow: ", nrow(res))
	message("ncpg: ", length(unique(res$outcome)))
	message("nexp: ", length(unique(res$exposure)))
	if(nrow(res) > 0)
	{
		res$nom <- nom[i]
		l[[i]] <- subset(res, pval < 1e-5, select=c(exposure, outcome, method, b, se, pval, nsnp))
	}
}

res <- bind_rows(l)
save(res,file="tophits_followup.rdata")

