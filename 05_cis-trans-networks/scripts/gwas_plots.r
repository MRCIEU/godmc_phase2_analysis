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

stab <- subset(dat, !grepl("metabolites__", fn) & nsnp > 2, select=c(clust, label, nsnp, min_p, fisher)) %>% arrange(fisher)
write.csv(stab, file="../results/gwas_clusters_nochr6.csv")


group_by(dat, clust) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% as.data.frame

dat <- subset(dat, nsnp != 0)
dat_sig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 < 0.05/nrow(dat))
dat_nsig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 >= 0.05/nrow(dat))

dat$fdr <- p.adjust(dat$fisher)

dat_sig <- subset(dat, !grepl("metabolites__", fn) & fisher < 0.05/nrow(dat) & nsnp > 2)
dat_nsig <- subset(dat, !grepl("metabolites__", fn) & fisher >= 0.05/nrow(dat) & nsnp > 2)
dat_sig$fisher[dat_sig$fisher < 1e-50] <- 1e-50

p1 <- ggplot(dat_nsig, aes(x=label, y=-log10(fisher))) +
geom_point(aes(size=nsnp)) +
geom_point(data=dat_sig, aes(colour=as.factor(clust), size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment", size="Number\nof SNPs in\ncommunity", colour="Cluster") +
scale_colour_brewer(type="qual") +
facet_grid(. ~ subcategory, scale="free", space="free") +
theme(legend.position=c(0.03, 0.76), strip.text=element_text(angle=90, size=10), axis.text.x=element_text(size = 10))
# geom_label_repel(data=dat_sig, aes(label=clust))
p1
ggsave(p1, file="../images/gwas_clusters_full.pdf", width=20, height=13)

# temp <- group_by(dat, id=id.y) %>%
# summarise(
# 	nsnp=mean(nsnp),
# 	nsig=sum(fisher < 0.05),
# 	fdr=sum(fdr < 0.05),
# 	bonferroni=sum(bonferroni3 < 0.05),
# 	minp=min(min_p, na.rm=TRUE)
# ) %>% filter(bonferroni > 0)

# ggplot(subset(dat, !grepl("metabolites__", fn) & nsnp >= 2 & id.y %in% temp$id), aes(x=label, y=-log10(fisher))) +
# geom_point(aes(colour=clust, size=nsnp)) +
# theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
# geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
# labs(x="", y="Enrichment -log10(p-value)", size="Number\nof SNPs in\ncommunity") +
# scale_colour_continuous(guide=FALSE) +
# coord_flip()
# ggsave("../images/gwas_clusters_filtered.pdf", width=12, height=10)


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



