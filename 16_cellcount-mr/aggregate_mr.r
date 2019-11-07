library(dplyr)

resdir <- "results"
fn <- list.files(file.path(resdir, "scratch")) %>% grep("mr", ., value=TRUE)

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

mrivw <- bind_rows(lmrivw) %>% as_tibble()
het <- bind_rows(lhet) %>% as_tibble()
wr <- bind_rows(lwr) %>% as_tibble()

save(mrivw, het, wr, file="results/cellcount_mr.rdata")
