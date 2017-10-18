library(tidyverse)
library(matrixStats)

load("../data/trans_clumped.rdata")
load("../data/annotations.rdata")

load_permutations <- function()
{
	load("../results/matrix/m0.rdata")
	real <- mat
	load("../results/matrix/m2.rdata")
	real[1:10,1:10]
	mat[1:10,1:10]

	res <- reshape2::melt(real)

	avail <- list.files("../results/matrix")
	avail <- avail[avail != "m0.rdata"]

	for(i in avail)
	{
		message(i)
		load(paste0("../results/matrix/", i))
		nom <- paste0("p", i)
		res[[nom]] <- reshape2::melt(mat)$value
	}
	return(res)
}

res <- load_permutations()

# There is one set of counts for the real mQTLs, and ~1000 for the permuted SNP-CpG pairs.
# There are 2504 annotations
# There are 2504^2 possible pairs of annotations
# 


# Only counts that are far away from the distribution are going to be significant after multiple testing, so can afford to get rid of any rows for which the target

## Real = 0; else perm

y <- 0


# round 1 : remove if replicated

vars <- res[,1:2]
res <- res[,-c(1:2)]
res <- as.matrix(res)

sp <- split(1:nrow(res), as.factor(vars[,2]))
for(i in 1:length(sp))
{
	message(i)
	mat <- res[sp[[i]], ]
	lab <- vars[sp[[i]], ]
	save(mat, lab, file=paste0("../data/matrix/m", i, ".rdata"))
}



ind <- sapply(1:nrow(res), function(x)
{
	! res[x,y+1] %in% res[x,-c(y+1)]
})

res2 <- res[ind,]

ran <- sapply(1:nrow(res2), function(x)
{
	rank(res2[x,])[y+1]
})

vars$value <- res[,1]
# vars$rank[ind] <- ran
vars$m <- rowMeans(res[,-1])
vars$s <- rowSds(res[,-1])


for(i in 1:ncol(res))
{
	message(i)
	p1 <- pnorm(res[,i], vars$m, vars$s, low=TRUE)
	p2 <- pnorm(res[,i], vars$m, vars$s, low=FALSE)
	vars[[paste0("p", i)]] <- pmin(p1, p2)
}




p1 <- pnorm(res[,2], vars$m, vars$s, low=TRUE)
p2 <- pnorm(res[,2], vars$m, vars$s, low=FALSE)
vars$p1 <- pmin(p1, p2)

p1 <- pnorm(res[,3], vars$m, vars$s, low=TRUE)
p2 <- pnorm(res[,3], vars$m, vars$s, low=FALSE)
vars$p2 <- pmin(p1, p2)

p1 <- pnorm(res[,4], vars$m, vars$s, low=TRUE)
p2 <- pnorm(res[,4], vars$m, vars$s, low=FALSE)
vars$p3 <- pmin(p1, p2)

p1 <- pnorm(res[,5], vars$m, vars$s, low=TRUE)
p2 <- pnorm(res[,5], vars$m, vars$s, low=FALSE)
vars$p4 <- pmin(p1, p2)

save(vars, file="../results/matrix_vars.rdata")

##

a <- rep(0, 774)
for(i in 1:774)
{
	a[i] <- sum(vars[,6+i] < 1e-9)
}

sort(a)

b <- vars[vars[,7] < 1e-9, ]



vars2 <- vars
vars2$value <- res[,2]

vars2$p1 <- pnorm(vars2$value, vars2$m, vars2$s, low=TRUE)
vars2$p2 <- pnorm(vars2$value, vars2$m, vars2$s, low=FALSE)
vars2$p <- pmin(vars2$p1, vars2$p2)
sum(vars2$p < 1e-8)
min(vars2$p)


a <- subset(vars, p < 1e-9)

annos <- unique(c(a$Var1, a$Var2))

vars$score <- -log10(vars$p)
vars$score[vars$value < vars$m] <- vars$score[vars$value < vars$m] * -1

m <- tidyr::spread(subset(vars, select=c(Var1, Var2, score)), Var1, score)
m <- as.matrix(m[,-1])
m[!is.finite(m)] <- 300
heatmap(m)

save(m, file="../results/mat.rdata")

d <- dist(m)


library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
heatmap.2(m[1:1000,1:1000], tracecol=NA, col=my_palette, labRow=FALSE, labCol = FALSE, key=FALSE)


save(res, file="../results/matrix_rank.rdata")


# stopCluster(cl)

sort(res[res$rank == ncol(res)-3,-c(1:2)][2,])
sort(res[res$rank == ncol(res)-3,-c(1:2,ncol(res))][2,])
hist(as.numeric(res[res$rank == ncol(res)-3,-c(1:2,ncol(res))][2,]))


ress3 <- as.matrix(res[res$rank == 1 | res$rank == ncol(res)-3,-c(1,2,ncol(res))])
dif <- sapply(1:nrow(ress3), function(x) {
	a <- max(ress3[x,-c(1)])
	b <- min(ress3[x,-c(1)])
	dp <- ress3[x,1] - a
	dn <- b - ress3[x,1]
	return(c(dp, dn))
})
dif <- t(dif)
head(dif)
res2 <- subset(res, rank == 1 | rank == ncol(res)-3)
res2$dif <- dif[,1]
res2$dif[dif[,1] < 0] <- dif[dif[,1] < 0, 2]

p <- function(n)
{
	library(ggplot2)
	a <- data.frame(val=as.numeric(sig[n, -c(1,2,ncol(sig), ncol(sig)-1)]), what=c("Real", rep("Perm", ncol(sig)-5)))
	p1 <- ggplot(a, aes(x=val)) +
	geom_histogram() +
	geom_vline(data=subset(a, what=="Real"), aes(xintercept=val))
	return(p1)
}


sig <- subset(res2, dif > 200)
save(sig, file="../results/matrix_sig.rdata")

load("../results/matrix_sig.rdata")

p(1)
p(2)
p(3)
p(4)
p(5)



normal_test

a <- as.numeric(sig[2,-c(1,2,3,ncol(sig), ncol(sig)-1)])
hist(a)


p1 <- function(x)
{
	a <- as.numeric(sig[x,-c(1,2,3,ncol(sig), ncol(sig)-1)])
	print(hist(a))
}

a1 <- function(x)
{
	a <- as.numeric(sig[x,-c(1,2,3,ncol(sig), ncol(sig)-1)])
	return(a)
}

p1(10000)

shapiro.test(a1(3))

shap <- apply(sig[,-c(1,2,3,ncol(sig),ncol(sig)-1)], 1, function(x)
{
	shapiro.test(x)$p
})

outlier <- apply(sig[,-c(1,2,ncol(sig),ncol(sig)-1)],1,function(x))
{
	m <- mean(x[-1])
	s <- sd(x[-1])
	p <- pnorm(x[1], m, s, low=FALSE)
	return(p)
}

m <- rowMeans(sig[,-c(1,2,3,ncol(sig),ncol(sig)-1)])
s <- apply(sig[,-c(1,2,3,ncol(sig),ncol(sig)-1)], 1, sd)
p <- pnorm(sig[,4], m, s, low=FALSE)
hist(p)