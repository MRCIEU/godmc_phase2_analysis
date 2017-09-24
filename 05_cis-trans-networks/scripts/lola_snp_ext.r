message("lola snp ext")
library(LOLA)
library(dplyr)

load("../data/lola/snp_granges.rdata")

tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")

message("global")
ext_global_snp <- runLOLA(community_snps, mqtl_snps, tfbsdb)

save(ext_global_snp, file="../results/ext_global_snp.rdata")

message("communities")
ext_communities_snp <- list()
for(i in 1:length(community_snps_separate))
{
	message(i)
	ext_communities_snp[[i]] <- runLOLA(community_snps_separate[[i]], community_snps, tfbsdb, cores=5)
	ext_communities_snp[[i]][, "userSet"] <- rep(names(community_snps_separate)[i], nrow(ext_communities_snp[[i]]))
}
ext_communities_snp <- bind_rows(ext_communities_snp)
save(ext_communities_snp, file="../results/ext_snp_communities.rdata")
