library(ggplot2)
library(coloc)
library(simulateGP)

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
	return(map)
}

load("../data/ld.rdata")

r <- ld[[1]][[1]]
map <- ld[[1]][[2]]
map1 <- sample_causal_variants(map, 3, 0.3)
o1 <- simulate_summary_data(r, map1, 25000)
map2 <- map1
map2$b <- map2$b * 0.2
o2 <- simulate_summary_data(r, map2, 25000)

plot(o2$bhat)

plot(o1$bhat)

run_coloc <- function(ss1, ss2, coverage, method)
{
	ss1 <- ss1[ss1$selected,]
	ss2 <- ss2[ss1$selected,]

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
	} else {
		ss1 <- subset(ss1, include)
		ss2 <- subset(ss2, include)
	}

	d1 <- list(
		pvalues = ss1$pval[selected]
	)
	d2 <- list(
		)
}





ncausal <- 3
rsq_trait1 <- 0.4
nid <- 25000
radius <- 500000
S <- 0

params <- expand.grid(
	n = 28000,
	region = 1,
	ncausal = c(1, 2, 3),
	coverage = c(1, 0.8, 0.6, 0.4, 0.2, 0.1),
	rsq_trait1 = c(0, 0.3, 0.1, 0.01),
	scale_trait2 = c(1, 0.5, 0.1),
	method = c("fill", "sparse"),
	nsim = c(1:100)
) %>% as_tibble

