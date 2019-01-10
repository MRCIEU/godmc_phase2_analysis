library(ggplot2)
library(dplyr)
library(ggrepel)
library(TwoSampleMR)
library(gridExtra)
library(tidyr)



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


dat_sig <- subset(dat, p < 0.05/nrow(dat))
dat_nsig <- subset(dat, p >= 0.05/nrow(dat))

p1 <- ggplot(dat %>% subset(!grepl("Difference", trait)), aes(x=trait, y=-log10(p))) +
geom_point(aes(size=nsnp.x)) +
geom_point(data=dat_sig, aes(colour=lor, size=nsnp.x)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="-log10(p) enrichment", size="Number\nof regions in\ncommunity", colour="log(OR)") +
# scale_colour_brewer(type="qual") +
facet_grid(. ~ subcategory, scale="free", space="free") +
theme(legend.position="none", strip.text=element_text(angle=90, size=10), axis.text.x=element_text(size = 8)) +
geom_label_repel(data=dat_sig, aes(label=clust), size=2)
p1
ggsave(p1, file="../images/gwas_clusters_full.pdf", width=18, height=13)
ggsave(p1, file="../images/gwas_clusters_full.png", width=18, height=13)


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

