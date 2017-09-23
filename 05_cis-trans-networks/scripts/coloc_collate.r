library(dplyr)

fn <- paste0("../results/coloc", 1:999, ".rdata")
l <- list()
for(i in fn)
{
	message(i)
	if(file.exists(i))
	{
		load(i)
		l[[i]] <- res
	} else {
		message(i, " not there")
	}
}
res <- bind_rows(l)

sum(res$H4 > 0.8)
save(res, file="../results/coloc.rdata")
