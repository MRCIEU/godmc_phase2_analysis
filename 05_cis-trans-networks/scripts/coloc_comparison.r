library(dplyr)

# Number of cis-trans pairs before coloc
load("../results/creg_tcpg.rdata")
length(unique(paste(dat$creg, dat$tcpg)))
# 1674157

# Number of cis-trans pairs after coloc
load("../results/graph_unpruned.rdata")
length(unique(paste(resdat$creg, resdat$tcpg)))
# 272608
