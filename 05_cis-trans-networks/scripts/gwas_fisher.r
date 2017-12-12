library(dplyr)
library(data.table)
library(TwoSampleMR)


load("../results/communities.rdata")

comms <- group_by(communities, cluster) %>% summarise(nsnp=length(unique(snp)), ncpg=length(unique(cpg))) %>% arrange(desc(nsnp)) %>% filter(nsnp > 5)

selcom <- subset(communities, cluster %in% comms$cluster) %>% arrange(cluster)

# convert to rsids


snp_1kg <- fread("../../10_mr-cpg-gwas/data/eur.bim.orig")

snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)

snp_1kg <- subset(snp_1kg, c1 == 1 & c2 == 1)
snp_1kg$snp <- paste0("chr", snp_1kg$V1, ":", snp_1kg$V4, ":SNP")
snp_1kg <- subset(snp_1kg, !duplicated(snp))

selcom <- merge(selcom, subset(snp_1kg, select=c(V2, snp)), by="snp")


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
}

save(dat, selcom, res, file="../results/gwas_clusters.rdata")
