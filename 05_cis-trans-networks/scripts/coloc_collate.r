library(dplyr)

fn <- paste0("../results/coloc", 1:414, ".rdata")
l <- list()
for(i in fn)
{
	message(i)
	load(i)
	l[[i]] <- res
}
res <- bind_rows(l)

sum(res$H4 > 0.8)
save(res, file="../results/coloc.rdata")
