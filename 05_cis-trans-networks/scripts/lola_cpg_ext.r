message("lola cpg ext")
library(LOLA)
library(dplyr)

load("../data/lola/cpg_granges.rdata")

tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")

message("global")
ext_global_cpg <- runLOLA(community_cpgs, mqtl_cpgs, tfbsdb)

message("creg")
ext_creg <- runLOLA(community_creg, mqtl_cpgs, tfbsdb)

message("tcpg")
ext_tcpg <- runLOLA(community_tcpg, mqtl_cpgs, tfbsdb)

save(ext_global_cpg, ext_creg, ext_tcpg, file="../results/ext_global_cpg.rdata")

message("communities")
ext_communities_cpg <- runLOLA(community_cpgs_separate, community_cpgs, tfbsdb, cores=5)
save(ext_communities_cpg, file="../results/ext_communities_cpg.rdata")

ext_communities_cpg$fdr <- p.adjust(10^-ext_communities_cpg$pValueLog, "fdr")

ext_communities_cpg_tophits <- group_by(ext_communities_cpg, userSet) %>%
	mutate(fdr2 = p.adjust(10^(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.05)
save(ext_communities_cpg_tophits, file="../results/ext_communities_cpg_tophits.rdata")

rm(ext_communities_cpg)
gc()

message("permutations")
ext_communities_cpg_perm <- runLOLA(community_cpgs_separate_perm, community_cpgs, tfbsdb, cores=5)
save(ext_communities_cpg_perm, file="../results/ext_communities_cpg_perm.rdata")


