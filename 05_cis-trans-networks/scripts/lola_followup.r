message("Running SNP enrichments")
library(dplyr)
library(LOLA)
library(igraph)
library(GenomicRanges)
library(htmltab)

load("../results/creg_tcpg.rdata")
load("../results/graph.rdata")
load("../results/communities.rdata")
load("../results/lola_global.rdata")
load("../results/lola_communities_tophits.rdata")
load("../results/lola_communities.rdata")
load("../../results/16/16_clumped.rdata")
load("../data/grinfo.rdata")
load("../data/snpcontrolsets_selection.rdata")
grinfo$snp <- paste(grinfo$snpchr, grinfo$snppos, grinfo$snptype, sep=":")

ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))

zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())

clumped <- subset(clumped, cpg %in% zhou & pval < 1e-8)
grinfo <- subset(grinfo, snp %in% clumped$snp)
grinfo <- inner_join(grinfo, ldinfo, by=c("snp"="SNP"))
grinfo2 <- subset(grinfo, !duplicated(snp))

temp1 <- inner_join(clumped, ldinfo, by=c("snp"="SNP"))

# DUE TO MEMORY ISSUES ONLY USING 50k mQTLs as background
index <- sample(1:nrow(temp1), 50000, replace=FALSE)
snpuniverse <- GRanges(seqnames=temp1$snpchr[index], ranges=IRanges(temp1$min[index], temp1$max[index]), strand="+")
community_universe <- GRanges(seqnames=grinfo2$snpchr, ranges=IRanges(grinfo2$min, grinfo2$max), strand="+")
names(community_universe) <- grinfo2$snp

temp <- lapply(split(enr_communities_tophits, enr_communities_tophits$antibody), function(x) {
	a <- subset(grinfo, cluster %in% x$userSet)
	a <- subset(a, !duplicated(snp))
	b <- GRanges(seqnames=a$snpchr, ranges=IRanges(a$min, a$max), strand="+")
	names(b) <- a$snp
	return(b)
}) %>% GRangesList


## Saving for future analysis
community_snps <- community_universe
community_snps_separate <- temp
save(community_snps, community_snps_separate, file="../data/lola/snp_granges.rdata")

tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

message("global")
enr_snp_global <- runLOLA(community_universe, snpuniverse, tfbsdb)

save(enr_snp_global, file="../results/lola_snp_global.rdata")

message("communities")
enr_snp_communities <- list()
for(i in 1:length(temp))
{
	message(i)
	 enr_snp_communities[[i]] <- runLOLA(temp[[i]], community_universe, tfbsdb, cores=5)
	 enr_snp_communities[[i]][, "userSet"] <- rep(names(temp)[i], nrow(enr_snp_communities[[i]]))
}
enr_snp_communities <- bind_rows(enr_snp_communities)
save(enr_snp_communities, file="../results/lola_snp_communities.rdata")

rm(tfbsdb)
gc()


tfbsdb2 <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")

message("global2")
enr_snp_global2 <- runLOLA(community_universe, snpuniverse, tfbsdb2)
save(enr_snp_global2, file="../results/lola_snp_global2.rdata")

message("communities2")
enr_snp_communities2 <- list()
for(i in 1:length(temp))
{
	message(i, " of ", length(temp))
	 enr_snp_communities2[[i]] <- runLOLA(temp[[i]], community_universe, tfbsdb2, cores=5)
	 enr_snp_communities2[[i]][, "userSet"] <- rep(names(temp)[i], nrow(enr_snp_communities2[[i]]))
}
enr_snp_communities2 <- bind_rows(enr_snp_communities2)
save(enr_snp_communities2, file="../results/lola_snp_communities2.rdata")



# Quick summary

# CpGs in communities tend to be at motifs related to cohesin - which relates to DNA organisation in the nucleosome.
# SNPs that influence CpG communities tend to be in regions to do with histone marks or RNA polymerase. 
# Pioneer transcription factors have a role in unwinding inaccessible DNA to initiate TF binding

# H3K9me3
# H3K4me1
# H3K36me3

# Use LD min/max for SNP region /panfs/panasas01/shared-godmc/1kg_reference_ph3/snpcontrolsets_selection.rdata
# For each CpG enrichment motif, what are the SNPs enriched for
# Add extended database which includes roadmap
# Look at differences across tissues for the same motif



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
