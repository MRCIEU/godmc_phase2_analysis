# 1. Filter communities so that there are no correlated CpGs due to physical proximity

message("lola cpg 25 states")

library(LOLA)
library(dplyr)
library(GenomicRanges)

load("../data/lola/cpg_granges.rdata")



library(tidyr)
a <- separate(s25_global_cpg, filename, c("celltype", "v1", "v2", "v3", "state"), sep="_")


# Read in stuff
regionDB <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

regionDB <- loadRegionDB("../../data/lola/",collection="25states")

message("global")
s25_global_cpg <- runLOLA(community_cpgs, mqtl_cpgs, regionDB)

message("creg")
s25_creg <- runLOLA(community_creg, mqtl_cpgs, regionDB)

message("tcpg")
s25_tcpg <- runLOLA(community_tcpg, mqtl_cpgs, regionDB)

save(s25_global_cpg, s25_creg, s25_tcpg, file="../results/s25_global_cpg.rdata")

message("communities")
s25_communities_cpg <- runLOLA(community_cpgs_separate, community_cpgs, regionDB, cores=5)
save(s25_communities_cpg, file="../results/s25_communities_cpg.rdata")

s25_communities_cpg$fdr <- p.adjust(10^-s25_communities_cpg$pValueLog, "fdr")

s25_communities_cpg_tophits <- group_by(s25_communities_cpg, userSet) %>%
	mutate(fdr2 = p.adjust(10^(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.05)
save(s25_communities_cpg_tophits, file="../results/s25_communities_cpg_tophits.rdata")

head(s25_communities_cpg_tophits)
dim(s25_communities_cpg_tophits)

rm(s25_communities_cpg)
gc()

message("permutations")
s25_communities_cpg_perm <- runLOLA(community_cpgs_separate_perm, community_cpgs, regionDB, cores=5)
save(s25_communities_cpg_perm, file="../results/s25_communities_cpg_perm.rdata")

rm(s25_communities_cpg_perm)
gc()

