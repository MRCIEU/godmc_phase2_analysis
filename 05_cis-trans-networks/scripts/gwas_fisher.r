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
	dat$nsig1[j] <- sum(pvals < 0.01)
	dat$nsig2[j] <- sum(pvals < 0.001)
	dat$nsig3[j] <- sum(pvals < 0.0001)
}

dat <- subset(dat, nsnp >= 5)

for(i in 1:nrow(dat))
{
	message(i)
	dat$binom1[i] <- binom.test(x=dat$nsig1[i], n=dat$nsnp[i], p=0.01)$p.v
	dat$binom2[i] <- binom.test(x=dat$nsig2[i], n=dat$nsnp[i], p=0.001)$p.v
	dat$binom3[i] <- binom.test(x=dat$nsig3[i], n=dat$nsnp[i], p=0.0001)$p.v
	dat$binom4[i] <- binom.test(x=dat$nsig3[i], n=dat$nsnp[i], p=0.001)$p.v
}


save(dat, selcom, res, file="../results/gwas_clusters.rdata")

library(ggplot2)
library(dplyr)
library(ggrepel)
load("../results/gwas_clusters_nochr6.rdata")

info <- read.csv("../../data/gwas/00info.csv")
info <- data.frame(fn=gsub(".txt.gz", "", info$newfile), id=info$id)

labels <- subset(res, !duplicated(i), select=c(i, fn))
labels$fn <- gsub("../../data/gwas/", "", labels$fn)
labels$fn <- gsub(".txt.gz", "", labels$fn)
dat <- merge(dat, labels, by.x="id", by.y="i")
dat$fisher <- unlist(dat$fisher)
dat$fdr <- p.adjust(dat$fisher, "fdr")
dat$bonferroni <- p.adjust(dat$fisher, "bonferroni")
dat$bonferroni3 <- p.adjust(dat$binom4, "bonferroni")
dat$label <- gsub("disease__", "", dat$fn)
dat$label <- gsub("risk_factor__", "", dat$label)
dat$label <- gsub("_", " ", dat$label)

stab <- subset(dat, select=c(clust, label, nsnp, min_p, nsig1, nsig2, nsig3, binom1, binom2, binom3, binom4, fisher)) %>% arrange(binom4)
write.csv(stab, file="../results/gwas_clusters_nochr6.csv")

dat <- merge(dat, info, by="fn")

library(TwoSampleMR)
ao <- available_outcomes()
#load("../data/outcomes.RData")
ao <- subset(ao, select=c(id, subcategory))
#dat <- merge(dat, ao, by.x="id.y", by.y="id")
dat <- merge(dat, ao, by.x="id.y", by.y="id",all.x=T)

dat$subcategory[dat$subcategory=="Hemodynamic"] <- "Haemotological"
dat$subcategory[dat$subcategory=="Immune system"] <- "Autoimmune / inflammatory"
dat$subcategory[dat$subcategory=="Diabetes"] <- "Glycemic"
dat$subcategory[dat$subcategory=="Biomarker"] <- "Other"
dat$subcategory[dat$subcategory=="Protein"] <- "Other"
dat$subcategory[dat$subcategory=="Reproductive aging"] <- "Aging"
dat$subcategory[dat$subcategory=="Lung disease"] <- "Other"
dat$subcategory[dat$subcategory=="Autoimmune / inflammatory"] <- "Immune"
dat$subcategory[dat$subcategory=="Psychiatric / neurological"] <- "Neurological"
dat$subcategory[is.na(dat$subcategory)] <- "Kidney"

group_by(dat, clust) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% as.data.frame


dat_sig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 < 0.05/nrow(dat))
dat_nsig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 >= 0.05/nrow(dat))

p1 <- ggplot(subset(dat_nsig, !grepl("metabolites__", fn) & nsnp > 3), aes(x=label, y=-log10(binom4))) +
geom_point(aes(size=nsnp)) +
geom_point(data=dat_sig, aes(colour=as.factor(clust), size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment", size="Number\nof SNPs in\ncommunity") +
scale_colour_brewer(type="qual", guide=FALSE) +
facet_grid(. ~ subcategory, scale="free", space="free") +
theme(legend.position="bottom", strip.text=element_text(angle=90, size=10), axis.text.x=element_text(size = 10)) +
geom_label_repel(data=dat_sig, aes(label=clust))
ggsave(p1, file="../images/gwas_clusters_full.pdf", width=20, height=13)

temp <- group_by(dat, id) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni3 < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% filter(bonferroni > 0)

ggplot(subset(dat, !grepl("metabolites__", fn) & nsnp >= 5 & id %in% temp$id), aes(x=label, y=-log10(binom4))) +
geom_point(aes(colour=clust, size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment -log10(p-value)", size="Number\nof SNPs in\ncommunity") +
scale_colour_continuous(guide=FALSE) +
coord_flip()
ggsave("../images/gwas_clusters_filtered.pdf", width=12, height=10)


sig <- subset(dat, bonferroni3 < 0.05)


subset(dat, label == "schizophrenia")

load("../results/graph.rdata")
load("../data/entrez_genes.rdata")

library(pathfindR)

table(mem$cpg %in% anno$ind)

i <- 8
message(i)
dir.create(paste0("../results/pathfindR/", i), recursive=TRUE)
setwd(paste0("../results/pathfindR/", i))

temp1 <- subset(anno, ind %in% subset(mem, cluster == i)$cpg) %>%
	filter(!duplicated(values), values != "")
temp1$ind <- 5
temp1$pval <- 1e-8
out[[j]] <- run_pathfindR(temp1, n_processes=1, pin_name = "IntAct")
out$cluster <- i
setwd("../../../../../scripts")


temp <- subset(selcom, cluster == 9 & !duplicated(V2))

table(temp$V2 %in% gtex_eqtl$SNP)



