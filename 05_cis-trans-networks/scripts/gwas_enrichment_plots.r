library(ggplot2)
library(dplyr)
library(ggrepel)
library(TwoSampleMR)
library(gridExtra)
library(tidyr)
library(qqman)
library(magrittr)



load("../results/gwas_enrichment.rdata")
load("../data/labelids.rdata")
ao <- available_outcomes()
#load("../data/outcomes.RData")
#dat <- merge(dat, ao, by.x="id.y", by.y="id")
dat <- merge(gwas_enrichment, ao, by="id")
dat <- subset(dat, access != "developer")
dat <- merge(dat, labelids, by="id")
dat$subcategory[dat$subcategory=="Hemodynamic"] <- "Haematological"
dat$subcategory[dat$subcategory=="Haemotological"] <- "Haematological"
dat$subcategory[dat$subcategory=="Immune system"] <- "Autoimmune / inflammatory"
dat$subcategory[dat$subcategory=="Diabetes"] <- "Glycemic"
dat$subcategory[dat$subcategory=="Biomarker"] <- "Other"
dat$subcategory[dat$subcategory=="Protein"] <- "Other"
dat$subcategory[dat$subcategory=="Hormone"] <- "Other"
dat$subcategory[dat$subcategory=="Reproductive aging"] <- "Aging"
dat$subcategory[dat$subcategory=="Lung disease"] <- "Other"
dat$subcategory[dat$subcategory=="Autoimmune / inflammatory"] <- "Immune"
dat$subcategory[dat$subcategory=="Psychiatric / neurological"] <- "Neurological"
dat$subcategory[is.na(dat$subcategory)] <- "Kidney"

dat$fdr <- p.adjust(dat$p, "fdr")

dat_sig <- subset(dat, fdr < 0.1)
dat_nsig <- subset(dat, fdr >= 0.1)

p1 <- ggplot(dat %>% subset(!grepl("Difference", trait)), aes(x=trait, y=-log10(p))) +
geom_point(aes(size=nsnp.x)) +
geom_point(data=dat_sig, aes(colour=lor, size=nsnp.x)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(max(dat_sig$p)), linetype="dotted") +
labs(x="", y="-log10(p) enrichment", size="Number\nof regions in\ncommunity", colour="log(OR)") +
# scale_colour_brewer(type="qual") +
facet_grid(. ~ subcategory, scale="free", space="free") +
theme(legend.position="none", strip.text=element_text(angle=90, size=10), axis.text.x=element_text(size = 8)) +
geom_label_repel(data=dat_sig, aes(label=clust), size=2)
p1
ggsave(p1, file="../images/gwas_clusters_full.pdf", width=18, height=13)
ggsave(p1, file="../images/gwas_clusters_full.png", width=18, height=13)

qq(dat$p)
median(qchisq(dat$p, 1,low=FALSE) / qchisq(0.5, 1))

o <- rep(NA, 1000)
for(i in 1:1000)
{
	o[i] <- median(qchisq(runif(nrow(dat)), 1,low=FALSE) / qchisq(0.5, 1))
}

quantile(o, 0.95)
hist(o)

####

pvals <- dat$p[is.finite(dat$p)]
temp <- data_frame(obs=-log10(sort(pvals)), exp=-log10((1:length(pvals) - 0.5)/(length(pvals))))
p2 <- ggplot(temp, aes(x=exp, y=obs)) +
geom_point() +
geom_abline(slope=1,intercept=0) +
labs(x="Expected -log10 (p)", y="Observed -log10 (p)")
p2


####

pp <- arrangeGrob(p2,p1,layout_matrix=
rbind(
	c(1,2,2,2),
	c(NA,2,2,2)
))
grid.arrange(pp)

ggsave("../images/gwas_enrichments.pdf", pp, width=18, height=10)
ggsave("../images/gwas_enrichments.png", pp, width=18, height=10)


####


temp <- spread(subset(dat, select=c(clust, trait, z)), key=clust, value=z)
rownames(temp) <- temp$trait
temp <- as.matrix(temp[,-1])
heatmap(temp)



####

# Are the regions for one community in preferential regions compared to other community regions?

load("../data/entity_info.rdata")
load("../data/snpcontrolsets_selection.rdata")

clusters <- unique(entities$cluster)
table(entities$snp_name %in%f.all$SNP)
f.all2 <- subset(f.all, SNP %in% entities$snp_name)
temp <- merge(f.all2, entities, by.x="SNP", by.y="snp_name")
l <- list()
for(i in 1:length(clusters))
{
	message(i)
	dum <- as.numeric(temp$cluster == clusters[i])
	o <- summary(glm(dum ~ nproxies+tssdist+GC_freq+CpG_freq+MAF, data=temp, family="binomial"))$coefficients[-1,]
	l[[i]] <- data_frame(cluster=clusters[i], n=sum(dum), factor=rownames(o), pval=o[,4])
}

enr_bias <- bind_rows(l)
save(enr_bias, file="../results/enr_bias.rdata")
load("../results/enr_bias.rdata")

ggplot(enr_bias, aes(x=pval)) +
geom_histogram() +
facet_grid(. ~ factor)
ggsave("../images/test_communit_enrichment_bias.pdf", width=10, height=5)

ggplot(subset(enr_bias, n >10), aes(x=pval)) +
geom_histogram() +
facet_grid(. ~ factor)
ggsave("../images/test_communit_enrichment_bias_gt10.pdf", width=10, height=5)

sig_clust <- subset(gwas_enrichment, p < 0.05/nrow(gwas_enrichment))$clust %>% unique

ggplot(subset(enr_bias, cluster %in% sig_clust), aes(x=pval)) +
geom_histogram() +
facet_grid(. ~ factor)
ggsave("../images/test_communit_enrichment_bias_gwassig.pdf", width=10, height=5)

temp <- subset(enr_bias, cluster %in% sig_clust)
subset(temp, p.adjust(pval, "fdr") < 0.05)

subset(enr_bias, cluster %in% sig_clust)$pval %>% p.adjust(., "fdr") %>% table(. < 0.05)


min(enr_bias$pval)
subset(enr_bias, pval < 1e-5)

enr_bias$fdr <- p.adjust(enr_bias$pval, "fdr")
subset(enr_bias,fdr < 0.05)



##

load("../results/core_communities_cpg_tophits.rdata")
load("../results/ext_communities_cpg_tophits.rdata")
load("../results/graph.rdata")
load("../data/entity_info.rdata")

subset(core_communities_cpg_tophits, userSet %in% dat_sig$clust)

subset(entities, cluster==32)

subset(core_communities_cpg_tophits, userSet %in% 32)
subset(ext_communities_cpg_tophits, userSet %in% 32) %>% as.data.frame

subset(ext_communities_cpg_tophits, userSet %in% dat_sig$clust) %>% group_by(collection) %>% summarise(n=n())


subset(core_communities_cpg_tophits, userSet %in% dat_sig$clust) %>% group_by(userSet, antibody) %>% summarise(n=n())

subset(gwas_enrichment, clust == 2) %$% min(p)