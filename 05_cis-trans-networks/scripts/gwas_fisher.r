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

library(ggplot2)
library(dplyr)
load("../results/gwas_clusters.rdata")

labels <- subset(res, !duplicated(i), select=c(i, fn))
labels$fn <- gsub("../../data/gwas/", "", labels$fn)
labels$fn <- gsub(".txt.gz", "", labels$fn)
dat <- merge(dat, labels, by.x="id", by.y="i")
dat$fisher <- unlist(dat$fisher)
dat$fdr <- p.adjust(dat$fisher, "fdr")
dat$bonferroni <- p.adjust(dat$fisher, "bonferroni")
dat$label <- gsub("disease__", "", dat$fn)
dat$label <- gsub("risk_factor__", "", dat$label)
dat$label <- gsub("_", " ", dat$label)

group_by(dat, clust) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% as.data.frame

ggplot(subset(dat, !grepl("metabolites__", fn) & nsnp > 3), aes(x=label, y=-log10(fisher))) +
geom_point(aes(colour=clust, size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment", size="Number\nof SNPs in\ncommunity") +
scale_colour_continuous(guide=FALSE)
ggsave("../images/gwas_clusters_full.pdf", width=20, height=13)

temp <- group_by(dat, id) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% filter(bonferroni > 0)

ggplot(subset(dat, !grepl("metabolites__", fn) & nsnp >= 5 & id %in% temp$id), aes(x=label, y=-log10(fisher))) +
geom_point(aes(colour=clust, size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment -log10(p-value)", size="Number\nof SNPs in\ncommunity") +
scale_colour_continuous(guide=FALSE)
ggsave("../images/gwas_clusters_filtered.pdf", width=12, height=10)
