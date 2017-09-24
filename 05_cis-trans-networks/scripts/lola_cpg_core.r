# 1. Filter communities so that there are no correlated CpGs due to physical proximity

message("lola cpg core")

library(LOLA)
library(dplyr)
library(GenomicRanges)

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
core_communities_cpg <- runLOLA(community_cpgs_separate, community_cpgs, tfbsdb, cores=5)
save(core_communities_cpg, file="../results/core_communities_cpg.rdata")

core_communities_cpg_tophits <- group_by(core_communities_cpg, userSet) %>%
	mutate(fdr2 = p.adjust(exp(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.05)
save(core_communities_cpg_tophits, file="../results/core_communities_cpg_tophits.rdata")

rm(core_communities_cpg)
gc()

message("permutations")
core_communities_cpg_perm <- runLOLA(community_cpgs_separate_perm, community_cpgs, tfbsdb, cores=5)
save(core_communities_cpg_perm, file="../results/core_communities_cpg_perm.rdata")

rm(core_communities_cpg_perm)
gc()

# Create files for SNP enrichments

load("../data/grinfo.rdata")
load("../data/lola/snp_granges.rdata")
community_snps_separate <- lapply(split(core_communities_cpg_tophits, core_communities_cpg_tophits$antibody), function(x) {
	a <- subset(grinfo, cluster %in% x$userSet)
	a <- subset(a, !duplicated(snp))
	b <- GRanges(seqnames=a$snpchr, ranges=IRanges(a$min, a$max), strand="+")
	names(b) <- a$snp
	return(b)
}) %>% GRangesList

save(mqtl_snps, community_snps, community_snps_separate, file="../data/lola/snp_granges.rdata")



###


# url <- "https://github.com/mdozmorov/genomerunner_web/wiki/ENCODE-cell-types"
# a <- htmltab(doc=url)


# ## ld dif
# proxytest <- rbind(subset(temp1, select=c(snp, min, max, nproxies)), subset(grinfo2, select=c(snp, min, max, nproxies)))
# proxytest$id <- 1
# proxytest$id[1:nrow(temp1)] <- 0

# summary(glm(id ~ nproxies, proxytest, family="binomial"))
# proxytest$dif <- proxytest$max - proxytest$min
# summary(glm(id ~ dif, proxytest, family="binomial"))
# tapply(proxytest$dif, proxytest$id, summary)
