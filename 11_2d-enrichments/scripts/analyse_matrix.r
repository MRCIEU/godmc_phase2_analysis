library(tidyverse)
load("../results/matrix/m0.rdata")

real <- mat
load("../results/matrix/m2.rdata")
real[1:10,1:10]
mat[1:10,1:10]



load("../data/trans_clumped.rdata")
load("../data/annotations.rdata")


res <- reshape2::melt(real)

for(i in 1:100)
{
	if(file.exists(paste0("../results/matrix/m", i, ".rdata")))
	{
		message(i)
		load(paste0("../results/matrix/m", i, ".rdata"))
		nom <- paste0("p", i)
		res[[nom]] <- reshape2::melt(mat)$value
	}
}

ord <- apply(res[,-c(1:2)], 1, function(x)
{
	n <- length(x)
	mid <- round(n/2)
	x1 <- x[c(2:mid, 1, (mid+1):n)]
	order(x1, decreasing=TRUE)[mid] / n
})

res$ord <- ord

sig <- subset(res, ord == 0.01)

pl <- function(x)

res <- data.frame(snp=res$Var1, cpg=res$Var2, count=res$value, emp_p=ord)
# Check columns are correct
# Try to get empirical pval based on how far away it is

