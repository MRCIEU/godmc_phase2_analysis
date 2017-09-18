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

snpuniverse <- GRanges(seqnames=temp1$snpchr, ranges=IRanges(temp1$min, temp1$max), strand="+")
community_universe <- GRanges(seqnames=grinfo2$snpchr, ranges=IRanges(grinfo2$min, grinfo2$max), strand="+")


tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")
tfbsdb2 <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLAExt/hg19")


enr_snp_global <- runLOLA(community_universe, snpuniverse, tfbsdb)y
enr_snp_global2 <- runLOLA(community_universe, snpuniverse, tfbsdb2)
save(enr_snp_global2, enr_snp_global, file="../results/lola_snp_global.rdata")

temp <- lapply(split(enr_communities_tophits, enr_communities_tophits$antibody), function(x) {
	a <- subset(grinfo, cluster %in% x$userSet)
	a <- subset(a, !duplicated(snp))
	return(GRanges(seqnames=a$snpchr, ranges=IRanges(a$min, a$max), strand="+"))
}) %>% GRangesList

enr_snp_communities <- runLOLA(temp, community_universe, tfbsdb, cores=10)
save(enr_snp_communities, file="../results/lola_snp_communities.rdata")

enr_snp_communities2 <- runLOLA(temp, community_universe, tfbsdb2, cores=10)
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
