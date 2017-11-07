library(tidyverse)
library(matrixStats)

load("../data/trans_clumped.rdata")
load("../data/annotations.rdata")

# There is one set of counts for the real mQTLs, and ~1000 for the permuted SNP-CpG pairs.
# There are 2504 annotations
# There are 2504^2 possible pairs of annotations
# 


# Only counts that are far away from the distribution are going to be significant after multiple testing, so can afford to get rid of any rows for which the target

## Real = 0; else perm

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

	for(i in avail[1:100])
	{
		message(i)
		load(paste0("../results/matrix/", i))
		nom <- paste0("p", i)
		res[[nom]] <- reshape2::melt(mat)$value
	}
	return(res)
}

get_ranks <- function(res, y)
{
	message("Subsetting vars")
	vars <- res[,1:2]
	res <- res[,-c(1:2)]
	res <- as.matrix(res)

	if(y != 0)
	{
		message("Testing for permutation ", y)
		res[,1] <- res[,y+1]
		res <- res[,-c(y+1)]
	} else {
		message("Running real")
	}

	message("Starting with ", nrow(res), " rows")

	# round 1 : remove row if real value is replicated
	message("Removing replicated")
	ind <- sapply(1:nrow(res), function(x)
	{
		! res[x,1] %in% res[x,-c(1)]
	})

	res2 <- res[ind,]
	message(nrow(res2), " rows remaining")


	# Get ranks amongst remaining
	message("Getting ranks")
	ran <- sapply(1:nrow(res2), function(x)
	{
		rank(res2[x,])[1]
	})

	# Get distance when rank is highest
	ress3 <- as.matrix(res2[ran == 1 | ran == ncol(res2),])

	message(nrow(ress3), " rows left with top rank")

	message("Getting distances")
	dif <- apply(ress3[,-1], 1, summary) %>% t %>% as.data.frame
	dif$val <- ress3[,1]

	dif$dif <- dif$val - dif$Max
	dif$dif[dif$dif < 0] <- dif$Min[dif$dif < 0] - dif$val[dif$dif < 0]
	dif$sds <- rowSds(ress3[,-1])
	dif$sddif <- abs(dif$val - dif$Mean) / dif$sds

	message("Removing invariants")
	# How often is everything invariant except real
	temp <- subset(dif, sds == 0)
	# What is the value of real and permutations when everything is invariant
	table(temp$val, temp$Mean)

	ind2 <- dif$sds == 0
	message(sum(ind2), " invariants")
	dif2 <- dif[!ind2, ]

	vars2 <- vars[ind,1:2]
	vars3 <- vars2[ran == 1 | ran == ncol(res2), ]
	vars4 <- vars3[!ind2,]

	difres <- cbind(vars4, dif2)
	return(difres)
}

args <- commandArgs(T)
y <- as.numeric(args[1])

res <- load_permutations()

difres <- get_ranks(res, y)

save(difres, file=paste0("../results/difres/difres", y, ".rdata"))







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






## Assume normal distribution

vars$value <- res[,1]
vars$m <- rowMeans(res[,-1])
vars$s <- rowSds(res[,-1])

sp <- split(1:nrow(res), as.factor(vars[,2]))
for(i in 1:length(sp))
{
	message(i)
	mat <- res[sp[[i]], ]
	lab <- vars[sp[[i]], ]
	save(mat, lab, file=paste0("../data/matrix/m", i, ".rdata"))
}

for(i in 1:ncol(res))
{
	message(i)
	p1 <- pnorm(res[,i], vars$m, vars$s, low=TRUE)
	p2 <- pnorm(res[,i], vars$m, vars$s, low=FALSE)
	vars[[paste0("p", i)]] <- pmin(p1, p2)
}

save(vars, file="../results/matrix_vars.rdata")




##

a <- rep(0, ncol(vars)-5)
for(i in 1:length(a))
{
	a[i] <- sum(vars[,5+i] < 1e-9)
}

sort(a)

b <- vars[vars[,6] < 1e-9, ]
annos <- unique(c(b$Var1, b$Var2))


vars$score <- -log10(vars$p1)
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