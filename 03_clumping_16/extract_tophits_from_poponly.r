library(dplyr)
library(parallel)
load("../results/16/16_clumped.rdata")

clumped <- ungroup(clumped)

o <- lapply(1:10, function(i)
{
	load("../results/16_poponly/16_1.")
})
