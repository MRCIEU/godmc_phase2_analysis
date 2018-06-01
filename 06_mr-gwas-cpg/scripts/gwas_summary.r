library(dplyr)

lambda <- function(x)
{
	x <- x[is.finite(x)]
	x[x == 0] <- 1e-300
	ntp <- length(x)
	x <- qchisq(x, 1, lower.tail=FALSE)
     x <- sort(x)
    ppoi <- ppoints(x)
    ppoi <- sort(qchisq(ppoi, df=1, lower.tail=FALSE))
    x <- x[1:ntp]
    ppoi <- ppoi[1:ntp]
    m <- median(x, na.rm=TRUE)/qchisq(0.5, df=1)
    s <- summary( lm(x~0+ppoi) )$coeff
    e <- s[1,1]
    se <- s[1,2]
    return(c(m,e,se))
}

load("../data/snps_gwas.rdata")
traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)

l <- list()
for(i in 1:length(traits))
{
	l[[i]] <- list()
}
for(i in 1:942)
{
	message(i)
	load(paste0("../results/out/out", i, ".rdata"))
	for(j in 1:length(traits))
	{
		cat(".")
		l[[j]][[i]] <- subset(res, exposure == traits[j])
	}
	cat("\n")
}

for(i in 1:length(traits))
{
	message(i)
	res <- bind_rows(l[[i]])
	save(res, file=paste0("../results/out/gwas", i, ".rdata"))
}


thresh <- 0.05 / 370000
l <- list()
for(i in 1:140)
{
	message(i)
	load(paste0("gwas", i, ".rdata"))
	m <- list()
	if(nrow(res) != 0)
	{
		m$trait <- res$exposure[1]
		m$ntest <- sum(!is.na(res$pval))
		m$nsig <- sum(res$pval < thresh, na.rm=T)
		m$nsigl <- sum(res$pval < 0.05, na.rm=T)
		o <- lambda(res$pval)
		m$lambda_med <- o[1]
		m$lambda_reg <- o[2]
		m$lambda_regse <- o[3]
		m$minpval <- min(res$pval, na.rm=T)		
	}
	l[[i]] <- as.data.frame(m, stringsAsFactors=FALSE)
}


s <- bind_rows(l)
save(s, file="../../results/gwas_summary.rdata")


a1 <- subset(a, exposure %in% subset(s, nsig != 0)$trait)
ids <- unique(a1$id.exposure)

temp <- list.files("../sig_dat")
temp2 <- sapply(temp, function(x) strsplit(x, split="\\.")[[1]][1])
table(ids %in% temp2)

l <- list()
for(i in 1:140)
{
	message(i)
	load(paste0("gwas", i, ".rdata"))
	if(nrow(res) != 0)
	{
		l[[i]] <- subset(res, pval < 1e-6)
	}
}
sig <- bind_rows(l)
save(sig, file="../../results/mrbase_sig.rdata")


