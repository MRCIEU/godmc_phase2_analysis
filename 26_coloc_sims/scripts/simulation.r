library(ggplot2)
library(coloc)
library(simulateGP)
library(parallel)

sample_causal_variants <- function(map, ncausal, h2, S=0, radius=500000)
{
	# get range given radius
	nsnp <- nrow(map)
	rang <- map$bp[nsnp] - map$bp[1]
	window <- c(map$bp[1] + radius, map$bp[nsnp] - radius)
	windowi <- c(which(map$bp > window[1])[1], which(map$bp > window[2])[1]-1)

	# get the first one somewhere within range
	causals <- rep(0, ncausal)
	causals[1] <- sample(windowi[1]:windowi[2], 1)

	# define window
	pos <- map$bp[causals[1]]
	window <- which(map$bp > (pos - radius) & map$bp < (pos + radius))

	# sample others from within window
	if(ncausal > 1)
	{
		causals[2:ncausal] <- sample(window[!window == causals[1]], ncausal - 1, replace=FALSE)
	}

	p <- simulateGP::generate_gwas_params(map$freq[causals], h2 = h2, S=S)
	map$selected <- 1:nrow(map) %in% window
	map$causal <- 1:nrow(map) %in% causals
	map$b <- 0
	map$b[map$causal] <- p$beta
	return(map)
}


simulate_summary_data <- function(r, map, nid)
{
	map$beta <- map$b %*% r %>% drop
	ss <- simulateGP::generate_gwas_ss(map$beta, map$freq, nid)
	map$bhat <- ss$bhat
	map$se <- ss$se
	map$pval <- ss$pval
	map$nid <- nid
	return(map)
}


run_coloc <- function(ss1, ss2, coverage, method)
{
	sel <- ss1$selected
	ss1 <- ss1[sel,]
	ss2 <- ss2[sel,]

	prob <- -log10(ss1$pval)
	prob[is.infinite(prob)] <- 300
	prob <- prob / max(prob)

	include <- sample(1:nrow(ss1), coverage*nrow(ss1), replace=FALSE, prob=prob)
	ss1$include <- ss2$include <- 1:nrow(ss1) %in% include

	if(method == "fill")
	{
		ss1$bhat[!ss1$include] <- 0
		ss1$se[!ss1$include] <- 1
		ss1$pval[!ss1$include] <- 1
	} else if(method == "sparse") {
		ss1 <- subset(ss1, include)
		ss2 <- subset(ss2, include)
	} else {
		stop("method: fill or sparse")
	}

	d1 <- list(
		pvalues = ss1$pval,
		N=ss1$nid,
		MAF=ss1$freq,
		beta=ss1$bhat,
		varbeta=ss1$se^2,
		sdY=1,
		type="quant"
	)
	d2 <- list(
		pvalues = ss2$pval,
		N=ss2$nid,
		MAF=ss2$freq,
		beta=ss2$bhat,
		varbeta=ss2$se^2,
		sdY=1,
		type="quant"
	)
	coloc.abf(d1, d2)
}

simulation <- function(param, ld)
{
	ld <- ld[[param$region]]
	map1 <- sample_causal_variants(ld$map, param$ncausal, param$rsq_trait1)
	ss1 <- simulate_summary_data(ld$ld, map1, param$nid)
	if(param$coloc)
	{
		map2 <- map1
		map2$b <- map2$b * sqrt(param$rsq_trait2) / sqrt(param$rsq_trait1)
	} else {
		map2 <- sample_causal_variants(ld$map, param$ncausal, param$rsq_trait2)
	}
	ss2 <- simulate_summary_data(ld$ld, map2, param$nid)
	param$coloc_result <- run_coloc(ss1, ss2, param$coverage, "sparse")$summary %>% {which.max(.[-1])}
	return(param)
}


load("../data/ld.rdata")

params <- expand.grid(
	nid = 20000,
	region = 1:2,
	ncausal = c(1, 3),
	coverage = c(1, 0.1, 0.05, 0.01, 0.005),
	rsq_trait1 = c(0.1, 0.01),
	rsq_trait2 = c(0.1, 0.01, 0),
	coloc = c(TRUE, FALSE),
	method = c("fill", "sparse"),
	nsim = c(1:500)
) %>% as_tibble

params

res <- mclapply(1:nrow(params), function(x) {
	message(x)
	tryCatch(
	{
		simulation(params[x,], ld)
	}, error = function(e) {
		params[x,]
	})
}, mc.cores=16) %>% bind_rows()

save(res, file="../data/coloc_results.rdata")
