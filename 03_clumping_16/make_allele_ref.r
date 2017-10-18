library(data.table)
bim <- fread("../data/ref/eur.bim")
bim <- subset(bim, select=c(V2, V5, V6))
names(bim) <- c("snp", "a1", "a2")
save(bim, file="../data/ref/alleles.rdata")

