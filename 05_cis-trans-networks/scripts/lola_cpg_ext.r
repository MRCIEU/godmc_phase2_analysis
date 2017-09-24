# 1. Filter communities so that there are no correlated CpGs due to physical proximity

message("lola cpg ext")

library(LOLA)

load("../data/lola/cpg_granges.rdata")

# Read in stuff
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

ext_communities_cpg_tophits <- group_by(ext_communities_cpg, userSet) %>%
	mutate(fdr2 = p.adjust(exp(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.5)
save(ext_communities_cpg_tophits, file="../results/ext_communities_cpg_tophits.rdata")

rm(ext_communities_cpg)
gc()

message("permutations")
ext_communities_cpg_perm <- runLOLA(community_cpgs_separate_perm, community_cpgs, tfbsdb, cores=5)
save(ext_communities_cpg_perm, file="../results/ext_communities_cpg_perm.rdata")


