library(dplyr)


l <- list()
for(i in 1:962)
{
	load(paste0("../results/16/16_", i, "_clumped.rdata"))
	l[[i]] <- clumped
}
clumped <- bind_rows(l)
names(clumped)[names(clumped) == "P-value"] <- "pval"

save(clumped, file="../results/16/16_clumped.rdata")

table(clumped$pval < 5e-8)
table(clumped$pval < 5e-14)

cs <- subset(clumped, pval < 5e-14)

table(cs$cis)
length(unique(cs$cpg))



