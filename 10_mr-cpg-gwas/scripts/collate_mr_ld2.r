library(dplyr)

fn <- read.csv("../../data/gwas/00info.csv")
param <- expand.grid(split=1:10, gwas=1:nrow(fn))
extnom <- paste0("../data/extracted/filtered_gwas_mqtl_conditional_ready_", fn$id[param$gwas], ".rdata")
param$id <- 1:nrow(param)

param$fn <- paste0("../results/mr_ld2/out_", param$id, ".rdata")
fnr <- list.files("../results/mr_ld2")
param <- subset(param, param$fn %in% paste0("../results/mr_ld2/",fnr))

dir.create("../results/mr_ld", show=F)
param$out <- paste0("../results/mr_ld/out_", param$gwas, ".rdata")
group_by(param, gwas) %>%
	do({
		x <- .
		message(x$gwas[1])
		l <- list()
		for(i in 1:nrow(x))
		{
			if(file.exists(x$fn[i]))
			{
				load(x$fn[i])
				l[[i]] <- res
			} else {
				message(i, " missing")
			}
		}
		res <- bind_rows(l)
		save(res, file=x$out[1])
		return(data.frame(n=nrow(res)))
	})

