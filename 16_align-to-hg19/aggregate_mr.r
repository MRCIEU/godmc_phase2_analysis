library(dplyr)

resdir <- "results"
fn <- list.files(file.path(resdir, "scratch"))

lmrivw <- list()
lhet <- list()
lwr <- list()
for(i in 1:length(fn))
{
	message(i)
	load(file.path(resdir, "scratch", fn[i]))
	lmrivw[[i]] <- mrivw
	lhet[[i]] <- het
	lwr[[i]] <- wr
}

mrivw <- bind_rows(lmrivw)
het <- bind_rows(lhet)
wr <- bind_rows(lwr)

save(mrivw, het, wr, file="results/cellcount_mr.rdata")
