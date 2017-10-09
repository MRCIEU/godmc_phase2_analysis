library(dplyr)
library(tidyr)
library(GenomicRanges)
library(igraph)

load("../../results/16/16_clumped.rdata")
load("../results/coloc.rdata")
load("../results/creg_tcpg.rdata")

zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())

resdat <- inner_join(dat, res, by=c("creg"="exposure", "tcpg"="outcome")) %>%
	filter(creg %in% zhou, tcpg %in% zhou) %>%
	filter(H4 > 0.8) %>%
	filter(creg_pval < 1e-10, tcpg_pval < 1e-14)
dim(resdat)

# gr <- graph_from_data_frame(resdat, directed=TRUE)
# plot(gr)
# wc <- cluster_walktrap(gr, steps=20)
# class(wc)
# length(wc)

# mem <- membership(wc)
save(resdat, file="../results/graph_unpruned.rdata")


## PRUNE THE GRAPH
## THERE ARE LOTS OF UNCLUMPED CIS EFFECTS

dat <- resdat
dat <- group_by(dat, creg) %>%
	mutate(creg_count = n())
dat <- group_by(dat, tcpg) %>%
	mutate(tcpg_count = n())
table(dat$tcpg_count)
table(dat$creg_count)


# Are there any creg-tcpg which are also tcpg-creg?
i1 <- paste(dat$creg, dat$tcpg)
i2 <- paste(dat$tcpg, dat$creg)
sum(i1 %in% i2)
sum(i2 %in% i1)

# Remove them
dim(dat)
dat <- subset(dat, ! i1 %in% i2)
dim(dat)

# For each tcpg-creg chromosome
# Get the most connected creg and its most connected tcpg = sentinel
# Remove every cpg pair with the creg and tcpg are both within 2.5mb of sentinel
# Repeat with the next top CpG
dat <- dplyr::group_by(dat, creg_chr, tcpg_chr) %>%
	arrange(desc(creg_count), desc(tcpg_count)) %>%
	do({
		x <- .
		if(nrow(x) == 1) return(x)
		l <- list()
		i <- 1
		while(nrow(x) > 0)
		{
			l[[i]] <- x[1,]
			i <- i + 1
			message(i, " : ", nrow(x))
			index <- (abs(x$creg_pos - x$creg_pos[1]) - 2500000) > 0 &
				(abs(x$tcpg_pos - x$tcpg_pos[1]) - 2500000) > 0
			x <- x[index, ]
		}
		x <- bind_rows(l)
		return(x)
	})




# Remove tcpg if it's 5mb within sentinal cpg
# dat <- dplyr::group_by(dat, creg, tcpg_chr) %>%
# 	arrange(desc(tcpg_count), tcpg_pval) %>%
# 	do({
# 		x <- .
# 		if(nrow(x) == 1) return(x)
# 		y <- x[(abs(x$tcpg_pos - x$tcpg_pos[1]) - 2500000) > 0, ]
# 		x <- rbind(x[1,], y)
# 		return(x)
# 		})

# clumped <- subset(clumped, cis & pval < 1e-10 & cpg %in% dat2$creg)
# dat3 <- group_by(dat2, tcpg, creg_chr) %>%
# 	arrange(desc(creg_count), creg_pval) %>%
# 	do({
# 		x <- .
# 		a <- subset(clumped, cpg %in% x$creg[1])
# 		x <- subset(x, snp %in% a$snp)
# 		# if(nrow(x) == 1) return(x)
# 		# y <- x[(abs(x$creg_pos - x$creg_pos[1]) - 2500000) > 0, ]
# 		# x <- rbind(x[1,], y)
# 		return(x)
# 		})


# dat <- group_by(dat, tcpg, creg_chr) %>%
# 	arrange(desc(creg_count), creg_pval) %>%
# 	do({
# 		x <- .
# 		if(nrow(x) == 1) return(x)
# 		y <- x[(abs(x$creg_pos - x$creg_pos[1]) - 1000000) > 0, ]
# 		x <- rbind(x[1,], y)
# 		return(x)
# 		})




# gr2 <- graph_from_data_frame(dat2, directed=TRUE)
# gr3 <- graph_from_data_frame(dat3, directed=TRUE)
gr <- graph_from_data_frame(dat, directed=TRUE)

wc <- cluster_walktrap(gr, steps=20)

# # Simplifying doesn't work:
# gr2a <- igraph::simplify(gr2)
# is_simple(gr2)
# is_simple(gr2a)

# gr2adf <- igraph::as_data_frame(gr2a)

# wc <- cluster_walktrap(gr, steps=20)
# class(wc)
# length(wc)

mem <- membership(wc)
mem <- data_frame(cpg=names(mem), cluster=as.numeric(mem))

save(dat, gr, wc, mem, file="../results/graph.rdata")


#######

temp <- bind_rows(
	data_frame(cpg=dat$creg, snp=dat$snp, type="creg"),
	data_frame(cpg=dat$tcpg, snp=dat$snp, type="tcpg")) %>% 
filter(!duplicated(paste(cpg, type)))
temp$id <- 1:nrow(temp)

temp <- inner_join(temp, mem, by="cpg")

# Check that position is the same as lola expects
load("../../data/lola/annotated_cpgs.RData")
temp <- merge(dat, cpgs_annot, by.x="creg", by.y="ID")
table(temp$creg_pos == temp$start)

rm(temp)
gc()



######

####### FORMAT FOR GRANGES

load("../results/graph.rdata")
grinfo <- data_frame(
	cpg=c(dat$creg, dat$tcpg),
	chr=c(dat$creg_chr, dat$tcpg_chr),
	pos=c(dat$creg_pos, dat$tcpg_pos),
	snp=c(dat$snp, dat$snp)
) %>% filter(!duplicated(cpg))
grinfo <- inner_join(grinfo, mem, by="cpg") %>% arrange(cluster)

grinfo <- tidyr::separate(grinfo, snp, c("snpchr", "snppos", "snptype"), remove=FALSE)
grinfo$snppos <- as.numeric(grinfo$snppos)

load("../data/snpcontrolsets_selection.rdata")
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
grinfo <- inner_join(grinfo, ldinfo, by=c("snp"="SNP"))

grinfo <- group_by(grinfo, cluster) %>%
	mutate(cluster_size = n())

save(grinfo, file="../data/grinfo.rdata")



##### CPG GRANGES


load("../../results/16/16_clumped.rdata")
clumped1 <- subset(clumped, pval < 1e-14)
clumped1 <- subset(clumped1, !duplicated(cpg))

universe1 <- data_frame(chr=clumped1$cpgchr, start=clumped1$cpgpos, end=clumped1$cpgpos, cpg=clumped1$cpg, pval=clumped1$pval, cis=clumped1$cis)
universe2 <- with(grinfo, data_frame(chr=chr, start=pos, end=pos, cpg=cpg, ID=cluster))
universe3 <- universe2
universe3$ID <- sample(universe3$ID)

community_cpgs_separate <- lapply(split(universe2, universe2$ID), function(x) {
	GRanges(seqnames=x$chr, ranges=IRanges(x$start, x$end), strand="*")
}) %>% GRangesList
community_cpgs_separate_perm <- lapply(split(universe3, universe3$ID), function(x) {
	GRanges(seqnames=x$chr, ranges=IRanges(x$start, x$end), strand="*")
}) %>% GRangesList


community_creg <- with(subset(dat, !duplicated(creg)), GRanges(seqnames=creg_chr, ranges=IRanges(creg_pos, creg_pos), strand="+"))
names(community_creg) <- subset(dat, !duplicated(creg))$creg
community_tcpg <- with(subset(dat, !duplicated(tcpg)), GRanges(seqnames=tcpg_chr, ranges=IRanges(tcpg_pos, tcpg_pos), strand="+"))
names(community_tcpg) <- subset(dat, !duplicated(tcpg))$tcpg


mqtl_cpgs <- GRanges(seqnames=universe1$chr, ranges=IRanges(universe1$start, universe1$end), strand="*")
names(mqtl_cpgs) <- universe1$cpg
community_cpgs <- GRanges(seqnames=universe2$chr, ranges=IRanges(universe2$start, universe2$end), strand="*")
names(community_cpgs) <- universe2$cpg

save(mqtl_cpgs, community_cpgs, community_cpgs_separate, community_cpgs_separate_perm, community_creg, community_tcpg, file="../data/lola/cpg_granges.rdata")


### SNP GRANGES


clumped1 <- subset(clumped, pval < 1e-8)
grinfo2 <- subset(grinfo, !duplicated(snp))

temp1 <- inner_join(clumped1, ldinfo, by=c("snp"="SNP"))

# DUE TO MEMORY ISSUES ONLY USING 50k mQTLs as background

set.seed(100)
index <- sample(1:nrow(temp1), 50000, replace=FALSE)
mqtl_snps <- GRanges(seqnames=temp1$snpchr[index], ranges=IRanges(temp1$min[index], temp1$max[index]), strand="+")
names(mqtl_snps) <- temp1$snp[index]
community_snps <- GRanges(seqnames=grinfo2$snpchr, ranges=IRanges(grinfo2$min, grinfo2$max), strand="+")
names(community_snps) <- grinfo2$snp


save(mqtl_snps, community_snps, file="../data/lola/snp_granges.rdata")


