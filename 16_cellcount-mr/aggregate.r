library(dplyr)
library(GenomicRanges)

resdir <- "results"
fn <- list.files(file.path(resdir, "scratch"))
res <- lapply(fn, function(i)
{
	message(i)
	load(file.path(resdir, "scratch", i))
	chrompos <- as_tibble(chrompos)
	chrompos$variant <- paste0(chrompos$seqnames, ":", chrompos$end - 250000)

	out <- lapply(names(l), function(n)
	{
		a <- lapply(l[[n]], function(r) 
		{
			if(class(r) == "try-error")
			{
				return(rep(NA, 6))
			} else {
				return(r)
			}
		})
		x <- do.call(rbind, a) %>% 
		as_tibble() %>%
		mutate(trait=n)
		return(cbind(chrompos, x))
	}) %>% bind_rows()
}) %>% bind_rows() %>% as_tibble()

save(res, file="results/cellcount_coloc.rdata")

