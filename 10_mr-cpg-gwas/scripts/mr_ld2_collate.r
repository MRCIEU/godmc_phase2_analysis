library(dplyr)

fn <- read.csv("../../data/gwas/00info.csv")
fn$row <- 1:nrow(fn)

## Collate the MR analyses

l <- list()
for(i in 1:nrow(fn))
{
	message(i)
	f <- paste0("../results/mr_ld/ext/temp/mr_ld/out_", i, ".rdata")
	if(file.exists(f))
	{
		load(f)
		if(nrow(res) > 0)
		{
			if("maxnsnp" %in% names(res))
			{
				res <- subset(res, select=-c(maxnsnp)) %>% as_data_frame
			} else {
				res <- ungroup(res)
			}
			l[[i]] <- subset(res, Pvalue < 1e-6)			
		} else {
			message("no rows")
		}
	} else {
		message("Missing")
	}
}

res <- bind_rows(l)
res <- inner_join(res, fn, by=c("outcome"="row"))
dim(res)

save(res, file="../results/mr_ld_tophits.rdata")

