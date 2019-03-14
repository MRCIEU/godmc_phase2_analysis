library(ggplot2)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- subset(ann, select=c(chr, pos))

load("../results/out/gwas32.rdata")
ann$cpg <- rownames(ann)
res <- merge(res, ann, by.x="outcome", by.y="cpg")
res$chr <- gsub("chr", "", res$chr)
res$chr[res$chr == "X"] <- 23
res$chr[res$chr == "Y"] <- 24
res$chr <- as.numeric(res$chr)

ggplot(res %>% as.data.frame, aes(x=pos, y=-log10(pval) * sign(b))) +
geom_point() +
geom_hline(yintercept=-log10(0.05/nrow(res))) +
geom_hline(yintercept=log10(0.05/nrow(res))) +
facet_grid(. ~ chr, space="free", scale="free") +
theme_bw() +
theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
labs(x="", y="Signed -log10 p-value of causal effect")
ggsave("../images/manhattan_menarche.png", width=15, height=6)

load("../results/out/gwas100.rdata")
ann$cpg <- rownames(ann)
res <- merge(res, ann, by.x="outcome", by.y="cpg")
res$chr <- gsub("chr", "", res$chr)
res$chr[res$chr == "X"] <- 23
res$chr[res$chr == "Y"] <- 24
res$chr <- as.numeric(res$chr)
res$pval[res$pval < 1e-40] <- 1e-40
ggplot(res %>% as.data.frame, aes(x=pos, y=-log10(pval) * sign(b))) +
geom_point() +
geom_hline(yintercept=-log10(0.05/nrow(res))) +
geom_hline(yintercept=log10(0.05/nrow(res))) +
facet_grid(. ~ chr, space="free", scale="free") +
theme_bw() +
theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
labs(x="", y="Signed -log10 p-value of causal effect")
ggsave("../images/manhattan_cpd.png", width=15, height=6)

load("../results/out/gwas102.rdata")
ann$cpg <- rownames(ann)
res <- merge(res, ann, by.x="outcome", by.y="cpg")
res$chr <- gsub("chr", "", res$chr)
res$chr[res$chr == "X"] <- 23
res$chr[res$chr == "Y"] <- 24
res$chr <- as.numeric(res$chr)
res$pval[res$pval < 1e-40] <- 1e-40
ggplot(res %>% as.data.frame, aes(x=pos, y=-log10(pval) * sign(b))) +
geom_point() +
geom_hline(yintercept=-log10(0.05/nrow(res))) +
geom_hline(yintercept=log10(0.05/nrow(res))) +
facet_grid(. ~ chr, space="free", scale="free") +
theme_bw() +
theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
labs(x="", y="Signed -log10 p-value of causal effect")
ggsave("../images/manhattan_faod.png", width=15, height=6)

library(TwoSampleMR)
ao <- available_outcomes()
exp <- subset(ao, grepl("Cigaret", trait))$id %>% extract_instruments()
subset(ao, grepl("menarche", trait))
out <- extract_outcome_data(exp$SNP, 1095)
dat <- harmonise_data(exp, out)

mr(dat)

exp <- extract_instruments(1095)
out <- extract_outcome_data(exp$SNP, 961)
dat <- harmonise_data(exp, out)
mr(dat)

