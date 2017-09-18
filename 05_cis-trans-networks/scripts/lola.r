# 1. Filter communities so that there are no correlated CpGs due to physical proximity

library(dplyr)
library(LOLA)
library(igraph)
library(GenomicRanges)

load("../results/creg_tcpg.rdata")
load("../results/graph.rdata")

zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())

dat <- subset(dat, creg %in% zhou & tcpg %in% zhou)

dat <- group_by(dat, creg) %>%
	mutate(creg_count = n())
dat <- group_by(dat, tcpg) %>%
	mutate(tcpg_count = n())
table(dat$tcpg_count)
table(dat$creg_count)

# Remove tcpg if it's 5mb within sentinal cpg
dat2 <- group_by(dat, creg, tcpg_chr) %>%
	arrange(desc(tcpg_count), tcpg_pval) %>%
	do({
		x <- .
		if(nrow(x) == 1) return(x)
		y <- x[(abs(x$tcpg_pos - x$tcpg_pos[1]) - 2500000) > 0, ]
		x <- rbind(x[1,], y)
		return(x)
		})

gr <- graph_from_data_frame(dat2, directed=TRUE)
wc <- cluster_walktrap(gr, steps=20)
class(wc)
length(wc)

mem <- membership(wc)
mem <- data_frame(cpg=names(mem), cluster=as.numeric(mem))

temp <- bind_rows(
	data_frame(cpg=dat2$creg, snp=dat2$snp, type="creg"),
	data_frame(cpg=dat2$tcpg, snp=dat2$snp, type="tcpg")) %>% 
filter(!duplicated(paste(cpg, type)))
temp$id <- 1:nrow(temp)

temp <- inner_join(temp, mem, by="cpg")




grinfo <- data_frame(
	cpg=c(dat2$creg, dat2$tcpg),
	chr=c(dat2$creg_chr, dat2$tcpg_chr),
	pos=c(dat2$creg_pos, dat2$tcpg_pos),
	snp=c(dat2$snp, dat2$snp)
) %>% filter(!duplicated(cpg))
grinfo <- inner_join(grinfo, mem, by="cpg") %>% arrange(cluster)

grinfo <- tidyr::separate(grinfo, snp, c("snpchr", "snppos", "snptype"))
grinfo$snppos <- as.numeric(grinfo$snppos)
save(grinfo, file="../data/grinfo.rdata")

# Check that position is the same as lola expects
load("../../data/lola/annotated_cpgs.RData")
temp <- merge(dat, cpgs_annot, by.x="creg", by.y="ID")
table(temp$creg_pos == temp$start)

rm(temp)
gc()

# 2. Create bed files

#- one for all CpGs in communities against every other mQTL
#- one for CpGs in each community against every other community CpG

load("../../results/16/16_clumped.rdata")
clumped <- subset(clumped, cpg %in% zhou & pval < 1e-14)

universe1 <- data_frame(chr=clumped$cpgchr, start=clumped$cpgpos, end=clumped$cpgpos, cpg=clumped$cpg, pval=clumped$pval, cis=clumped$cis)
write.table(universe1, file="../data/lola/allmqtl.bed", row=F, col=F, qu=F, sep="\t")

universe2 <- with(grinfo, data_frame(chr=chr, start=pos, end=pos, cpg=cpg, ID=cluster))
write.table(universe2, file="../data/lola/communitymqtl.bed", row=F, col=F, qu=F, sep="\t")

universe3 <- universe2
universe3$ID <- sample(universe3$ID)
write.table(universe2, file="../data/lola/communitymqtl-permuted.bed", row=F, col=F, qu=F, sep="\t")

# group_by(universe2, ID) %>%
# do({
# 	write.table(., file=paste0("../data/lola/cluster", .$ID[1], ".bed"), row=F, col=F, qu=F)
# })

userset <- lapply(split(universe2, universe2$ID), function(x) {
	GRanges(seqnames=x$chr, ranges=IRanges(x$start, x$end), strand="+")
}) %>% GRangesList
userset_perm <- lapply(split(universe3, universe3$ID), function(x) {
	GRanges(seqnames=x$chr, ranges=IRanges(x$start, x$end), strand="+")
}) %>% GRangesList


allmqtl <- readBed("../data/lola/allmqtl.bed")
communitymqtl <- readBed("../data/lola/communitymqtl.bed")
cregmqtl <- with(subset(dat2, !duplicated(creg)), GRanges(seqnames=creg_chr, ranges=IRanges(creg_pos, creg_pos), strand="+"))
tcpgmqtl <- with(subset(dat2, !duplicated(tcpg)), GRanges(seqnames=tcpg_chr, ranges=IRanges(tcpg_pos, tcpg_pos), strand="+"))


tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

enr_global <- runLOLA(communitymqtl, allmqtl, tfbsdb)
enr_creg <- runLOLA(cregmqtl, allmqtl, tfbsdb)
enr_tcpg <- runLOLA(tcpgmqtl, allmqtl, tfbsdb)
save(enr_global, enr_creg, enr_tcpg, file="../results/lola_global.rdata")

enr_communities <- runLOLA(userset, communitymqtl, tfbsdb, cores=10)
save(enr_communities, file="../results/lola_communities.rdata")

enr_communities_tophits <- group_by(enr_communities, userSet) %>%
	mutate(fdr2 = p.adjust(exp(-pValueLog), "fdr")) %>%
	filter(fdr2 < 0.5)
save(enr_communities_tophits, file="../results/lola_communities_tophits.rdata")

enr_communities_perm <- runLOLA(userset_perm, communitymqtl, tfbsdb, cores=10)
save(enr_communities_perm, file="../results/lola_communities_perm.rdata")

#####

x <- subset(enr_communities, collection == "encode_tfbs")
length(unique(x$description))
length(unique(x$filename))
length(unique(x$userSet))
thresh <- -log10(0.05/nrow(x))
x <- subset(x, pValueLog > thresh)
table(x$description)

length(unique(x$userSet))

x1 <- subset(enr_communities, collection=="encode_tfbs" & pValueLog > -log10(0.05/length(unique(x$filename))))
length(unique(x1$userSet))
unique(x1$userSet)

filter(mem, cluster %in% unique(x1$userSet)) %>% group_by(cluster) %>% summarise(n=n()) %>% arrange(n) %>% as.data.frame


