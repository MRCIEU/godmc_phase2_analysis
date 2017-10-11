library(dplyr)

fn <- paste0("../results/annotations/annot_perm", 1:57, ".rdata")
l <- list()
for(i in fn)
{
	message(i)
	load(i)
	l[[i]] <- perms
}
perms <- bind_rows(l)

perms$real <- perms$perm == 0

out <- group_by(perms, mem) %>%
	do({
		x <- .
		sig <- sum(x$sig[!x$real] > x$sig[x$real]) / (nrow(x)-1)
		fdr <- sum(x$fdr[!x$real] > x$fdr[x$real]) / (nrow(x)-1)
		return(data.frame(sig=sig,fdr=fdr))
	})

save(perms, file="../results/annot_perms.rdata")




fn <- paste0("../results/annotations/annot", 1:57, ".rdata")
l <- list()
for(i in fn)
{
	message(i)
	load(i)
	l[[i]] <- subset(res, pval < 0.001)
}
res <- bind_rows(l)
dim(res)
save(res, file="../results/annot_tophits.rdata")
