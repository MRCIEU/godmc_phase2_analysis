library(dplyr)
library(data.table)
library(TwoSampleMR)

load("../data/cluster_extract.rdata")

res <- bind_rows(l)
dim(res)
table(res$i)

dat <- expand.grid(
	clust=unique(selcom$cluster), 
	id=1:length(unique(res$i)),
	nsnp=NA,
	min_p=NA,
	fisher=NA
)

for(j in 1:nrow(dat))
{
	message(j, " of ", nrow(dat))
	snps <- unique(subset(selcom, cluster == dat$clust[j])$V2)
	pvals <- na.omit(subset(res, i %in% dat$id[j] & V1 %in% snps)$V7)

	dat$nsnp[j] <- length(pvals)
	dat$min_p[j] <- min(pvals)
	dat$fisher[j] <- fishers_combined_test(pvals)$pval
	dat$nsig1[j] <- sum(pvals < 0.01)
	dat$nsig2[j] <- sum(pvals < 0.001)
	dat$nsig3[j] <- sum(pvals < 0.0001)
}

dat2 <- subset(dat, nsnp >= 4)

for(i in 1:nrow(dat2))
{
	message(i)
	dat2$binom1[i] <- binom.test(x=dat2$nsig1[i], n=dat2$nsnp[i], p=0.01)$p.v
	dat2$binom2[i] <- binom.test(x=dat2$nsig2[i], n=dat2$nsnp[i], p=0.001)$p.v
	dat2$binom3[i] <- binom.test(x=dat2$nsig3[i], n=dat2$nsnp[i], p=0.0001)$p.v
	dat2$binom4[i] <- binom.test(x=dat2$nsig3[i], n=dat2$nsnp[i], p=0.001)$p.v
}

dat <- merge(dat, subset(dat2, select=c(clust, id, binom1, binom2, binom3, binom4)), by=c("clust", "id"), all=TRUE)

save(dat, selcom, res, file="../results/gwas_clusters.rdata")

# Now without chr6

dat <- expand.grid(
	clust=unique(selcom$cluster), 
	id=1:length(unique(res$i)),
	nsnp=NA,
	min_p=NA,
	fisher=NA
)
selcom <- subset(selcom, !grepl("chr6", snp))

temp <- group_by(selcom, cluster) %>% summarise(n=length(unique(snp)))

for(j in 1:nrow(dat))
{
	message(j, " of ", nrow(dat))
	snps <- unique(subset(selcom, cluster == dat$clust[j])$V2)
	pvals <- na.omit(subset(res, i %in% dat$id[j] & V1 %in% snps)$V7)

	dat$nsnp[j] <- length(pvals)
	dat$min_p[j] <- min(pvals)
	dat$fisher[j] <- fishers_combined_test(pvals)$pval
	dat$nsig1[j] <- sum(pvals < 0.01)
	dat$nsig2[j] <- sum(pvals < 0.001)
	dat$nsig3[j] <- sum(pvals < 0.0001)
}

dat2 <- subset(dat, nsnp >= 4)

for(i in 1:nrow(dat2))
{
	message(i)
	dat2$binom1[i] <- binom.test(x=dat2$nsig1[i], n=dat2$nsnp[i], p=0.01)$p.v
	dat2$binom2[i] <- binom.test(x=dat2$nsig2[i], n=dat2$nsnp[i], p=0.001)$p.v
	dat2$binom3[i] <- binom.test(x=dat2$nsig3[i], n=dat2$nsnp[i], p=0.0001)$p.v
	dat2$binom4[i] <- binom.test(x=dat2$nsig3[i], n=dat2$nsnp[i], p=0.001)$p.v
}

dat <- merge(dat, subset(dat2, select=c(clust, id, binom1, binom2, binom3, binom4)), by=c("clust", "id"), all=TRUE)
save(dat, selcom, res, file="../results/gwas_clusters_nochr6.rdata")


