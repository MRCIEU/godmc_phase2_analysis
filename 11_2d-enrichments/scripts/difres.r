# library(tidyverse)
library(dplyr)
library(matrixStats)

# load("../data/trans_clumped.rdata")
# load("../data/annotations.rdata")

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


