library(dplyr)
a <- readRDS("../data/hm450.hg19.manifest.rds")
b <- as.data.frame(a)
b$cpg <- rownames(b)
load("../results/graph.rdata")

creg <- unique(dat$creg)
b$creg <- b$cpg %in% creg

table(b$creg, b$MASK_sub25_copy) %>% fisher.test

