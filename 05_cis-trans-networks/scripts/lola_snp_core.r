message("lola snp core")
library(LOLA)


load("../data/lola/snp_granges.rdata")

tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

message("global")
core_global_snp <- runLOLA(community_snps, mqtl_snps, tfbsdb)

save(core_global_snp, file="../results/core_global_snp.rdata")

message("communities")
core_communities_snp <- list()
for(i in 1:length(community_snps_separate))
{
	message(i)
	core_communities_snp[[i]] <- runLOLA(community_snps_separate[[i]], community_snps, tfbsdb, cores=5)
	core_communities_snp[[i]][, "userSet"] <- rep(names(community_snps_separate)[i], nrow(core_communities_snp[[i]]))
}
core_communities_snp <- bind_rows(core_communities_snp)
save(core_communities_snp, file="../results/core_snp_communities.rdata")
