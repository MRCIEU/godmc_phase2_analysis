library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

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

load("../data/traitlist.rdata")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- subset(ann, select=c(chr, pos))
chr6_cpgs <- rownames(ann)[ann$chr == "chr6"]
tr <- data.frame(id=1:length(traits), trait=traits)
lam <- expand.grid(chr6=c(TRUE, FALSE), id=1:length(traits), ncpg=NA, nsnp=NA, minp=NA, nsig=NA, lambda=NA, proppos=NA, pdir=NA)
lam <- merge(lam, tr, by="id")
for(i in 1:nrow(lam))
{
	message(i)
	if(lam$chr6[i])
	{
		load(paste0("../results/out/gwas", lam$trait[i], ".rdata"))
	} else {
		res <- subset(res, !outcome %in% chr6_cpgs)
	}
	if(nrow(res) > 0)
	{
		res$fdr <- p.adjust(res$pval, "fdr")
		lam$ncpg[i] <- nrow(res)
		lam$nsnp[i] <- max(res$nsnp)
		lam$minp[i] <- min(res$pval)
		lam$nsig[i] <- sum(res$fdr < 0.05)
		lam$lambda[i] <- lambda(res$pval)[1]
		lam$proppos[i] <- sum(sign(res$b)==1) / nrow(res)
		lam$pdir[i] <- binom.test(x=sum(sign(res$b) == 1), n=nrow(res), p=0.5)$p.value
	}
}
lam <- merge(lam, tr, by)
lam$phen <- do.call(c, lapply(as.character(lam$trait), function(x) strsplit(x, split="\\|")[[1]][1])) %>% gsub(" $", "", .)

lama <- subset(lam, chr6)
lamb <- subset(lam , !chr6)
lamb <- lamb %>% arrange(desc(lambda)) %>% subset(!duplicated(phen))
lamb$phen <- as.factor(lamb$phen)
lamb$phen <- factor(lamb$phen, levels=lamb$phen[order(lamb$lambda)])

lambtop <- subset(lamb, lambda > 1.05)

l <- list()
for(i in 1:nrow(lambtop))
{
	message(i)
	load(paste0("../results/out/gwas", lambtop$id[i], ".rdata"))
	res <- res[sample(1:nrow(res)), ]
	res$cut <- cut(1:nrow(res), 500)
	out <- res %>% group_by(cut) %>%
	summarise(lambda=lambda(pval)[1])
	out$id <- lambtop$id[i]
	out$phen <- lambtop$phen[i]
	l[[i]] <- out
}
lambdatop <- bind_rows(l)

ann$cpg <- rownames(ann)
l <- list()
for(i in 1:nrow(lambtop))
{
	message(i)
	load(paste0("../results/out/gwas", lambtop$id[i], ".rdata"))
	res <- merge(res, ann, by.x="outcome", by.y="cpg")
	out <- res %>% as_data_frame %>% group_by(chr) %>%
	summarise(n=n(), lambda=lambda(pval)[1])
	out$id <- lambtop$id[i]
	out$phen <- lambtop$phen[i]
	l[[i]] <- out
}
lambdatopchr <- bind_rows(l)


l <- list()
for(i in 1:nrow(lambtop))
{
	message(i)
	load(paste0("../results/out/gwas", lambtop$id[i], ".rdata"))
	res <- merge(res, ann, by.x="outcome", by.y="cpg")
	res <- subset(res, chr %in% paste("chr", 1:22, sep="")) %>% as_data_frame
	m <- list()
	for(j in 1:22)
	{
		temp <- subset(res, chr != paste0("chr", j))
		m[[j]] <- data_frame(chr=j, n=nrow(temp), lambda=lambda(temp$pval)[1])
	}
	out <- bind_rows(m)
	out$id <- lambtop$id[i]
	out$phen <- lambtop$phen[i]
	l[[i]] <- out
}
lambdatoploo <- bind_rows(l)

save(lam, lama, lamb, lambtop, lambdatop, lambdatopchr, lambdatoploo, file="../results/lambda.rdata")
