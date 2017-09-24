# 1. Filter communities so that there are no correlated CpGs due to physical proximity

library(LOLA)

load("../data/lola/cpg_granges.rdata")

# Read in stuff
tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

message("global")
core_global_cpg <- runLOLA(community_cpgs, mqtl_cpgs, tfbsdb)

message("creg")
core_creg <- runLOLA(community_creg, mqtl_cpgs, tfbsdb)

message("tcpg")
core_tcpg <- runLOLA(community_tcpg, mqtl_cpgs, tfbsdb)

save(core_global_cpg, core_creg, core_tcpg, file="../results/core_global_cpg.rdata")

message("communities")
core_communities <- runLOLA(userset, community_cpgs, tfbsdb, cores=5)
save(core_communities, file="../results/core_communities.rdata")

core_communities_tophits <- group_by(core_communities, userSet) %>%
	mutate(fdr2 = p.adjust(exp(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.5)
save(core_communities_tophits, file="../results/core_communities_tophits.rdata")

rm(core_communities)
gc()

message("permutations")
core_communities_perm <- runLOLA(userset_perm, community_cpgs, tfbsdb, cores=5)
save(core_communities_perm, file="../results/core_communities_perm.rdata")


