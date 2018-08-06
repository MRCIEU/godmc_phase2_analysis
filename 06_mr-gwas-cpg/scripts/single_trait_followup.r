
# BMI, HDL, LDL, trigs
# 

library(TwoSampleMR)
ao <- available_outcomes()
library(dplyr)

ewas <- read.table("../data/EWAS_Catalog_20-02-2018.txt.gz", he=T, sep="\t", stringsAsFactors=FALSE)
ewas$code <- paste(ewas$PMID, ewas$Trait)
b <- group_by(subset(ewas, P < 1e-7), code) %>% summarise(n=n()) %>% arrange(desc(n))

traits2 <- read.csv("../../data/gwas/00info.csv")
traits2$id_10 <- 1:nrow(traits2)
load("../data/snps_gwas.rdata")

ao$name <- paste(ao$trait, "||", ao$consortium, "||", ao$year, "||", ao$unit)
traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)
table(traits %in% ao$name)
rename <- traits[which(!traits %in% ao$name)]

a$exposure[a$exposure == rename[1]] <- "Serum creatinine (eGFRcrea) || CKDGen || 2015 || log ml/min/1.73 m^2"
a$exposure[a$exposure == rename[2]] <- "Serum cystatin C (eGFRcys) || CKDGen || 2015 || log ml/min/1.73 m^2"
a$exposure[a$exposure == rename[3]] <- "Father's age at death || UK Biobank || 2016 || SD"
a$exposure[a$exposure == rename[4]] <- "Mother's age at death || UK Biobank || 2016 || SD"


a0 <- filter(a, data_source.exposure == "mrbase") %>% group_by(id.exposure, exposure) %>% summarise(n=n()) %>% arrange(exposure, n)
a0$id_06 <- 1:nrow(a0)

remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)

a0 <- subset(a0, !id_06 %in% remove_ids)

table(a0$exposure %in% ao$name)

temp <- merge(subset(traits2, select=c(id, id_10)), ao, by="id")
table(a0$exposure %in% temp$name)

a0 <- merge(subset(a0, select=c(exposure, id_06)), temp, by.x="exposure", by.y="name", all=TRUE)

temp <- subset(a0, !is.na(id_06) | !is.na(id_10))

# Some traits are duplicated, look by eye to remove the out of date ones
remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)
a1 <- a0[-remove_ids,]
data.frame(substr(as.character(a1$exposure), start=1, stop=50), a1$n)



# mrbase ID
# 06 ID
# 10 ID
trait_ids <- data_frame()




lambda <- function(x)
{
	x <- x[is.finite(x)]
	x[x == 0] <- 1e-300
	ntp <- length(x)
	x <- qchisq(x, 1, lower.tail=FALSE)
	x <- sort(x)
	ppoi <- ppoints(x)
	ppoi <- sort(qchisq(ppoi, df=1, lower.tail=FALSE))
	x <- x[1:ntp]
	ppoi <- ppoi[1:ntp]
	m <- median(x, na.rm=TRUE)/qchisq(0.5, df=1)
	s <- summary( lm(x~0+ppoi) )$coeff
	e <- s[1,1]
	se <- s[1,2]
	return(c(m,e,se))
}


load("../../10_mr-cpg-gwas/results/mr_ld/out_41.rdata")
lambda(res$Pvalue)



args <- commandArgs(T)
jid <- as.numeric(args[1])
load("../data/snps_gwas.rdata")
load("../results/gwas_summary.rdata")
traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)
traits2 <- read.csv("../../data/gwas/00info.csv")
traits2$id2 <- 1:nrow(traits2)

grep("Body mass", traits)
jid <- 7
cpglist <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/bmi_cpg.txt", what=character())
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)

min(res$pval)

grep("Cigaret", traits)
jid <- 100
load("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/joehanes.rdata")

load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)


x1 <- data_frame(cpg=joehanes$Probe.ID, eff=joehanes$Effect)
x2 <- data_frame(cpg=res$outcome, b=res$b)
x <- merge(x1, x2)

cor(x$eff, x$b)
summary(lm(eff ~ b, x))



grep("choles", traits)
jid <- 22
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
cpglist <- subset(ewas, code %in% "28194238 High-density lipoprotein cholesterol")$CpG %>% unique
lambda(subset(res, outcome %in% cpglist)$pval)
fishers_combined_test(subset(res, outcome %in% cpglist)$pval)


# Only 4
jid <- 23
load(paste0("../results/out/gwas", jid, ".rdata"))
cpglist <- subset(ewas, code %in% "28194238 Total cholesterol")$CpG %>% unique
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)
fishers_combined_test(subset(res, outcome %in% cpglist)$pval)


jid <- 24
load(paste0("../results/out/gwas", jid, ".rdata"))
cpglist <- subset(ewas, code %in% "28213390 Serum low-density lipoprotein cholesterol")$CpG %>% unique
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)
fishers_combined_test(subset(res, outcome %in% cpglist)$pval)


jid <- 21
load(paste0("../results/out/gwas", jid, ".rdata"))
cpglist <- subset(ewas, code %in% "28194238 Triglycerides")$CpG %>% unique
lambda(res$pval)
lambda(subset(res, outcome %in% cpglist)$pval)
fishers_combined_test(subset(res, outcome %in% cpglist)$pval)


jid <- 7
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
cpglist <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/bmi_cpg.txt", what=character())
lambda(subset(res, outcome %in% cpglist)$pval)
cpglist <- subset(ewas, code %in% "28002404 Body mass index")$CpG %>% unique
lambda(subset(res, outcome %in% cpglist)$pval)


grep("Cigaret", traits)
jid <- 100
load(paste0("../results/out/gwas", jid, ".rdata"))
lambda(res$pval)
smok <- read.csv("../data/smok.csv")
cpglist <- subset(smok, Dose.Response == 1)$X...Name
lambda(subset(res, outcome %in% cpglist)$pval)
load("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/joehanes.rdata")


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- subset(ann, select=c(chr, pos))
chr6_cpgs <- rownames(ann)[ann$chr == "chr6"]
tr <- data.frame(id=1:length(traits), trait=traits)
lam <- expand.grid(chr6=c(TRUE, FALSE), id=1:length(traits), ncpg=NA, nsnp=NA, minp=NA, nsig=NA, lambda=NA, proppos=NA, pdir=NA)
lam <- merge(lam, tr, by="id")
for(i in 1:nrow(lam))
{
	message(i)
	if(lam$chr6[i])
	{
		load(paste0("../results/out/gwas", lam$trait[i], ".rdata"))
	} else {
		res <- subset(res, !outcome %in% chr6_cpgs)
	}
	if(nrow(res) > 0)
	{
		res$fdr <- p.adjust(res$pval, "fdr")
		lam$ncpg[i] <- nrow(res)
		lam$nsnp[i] <- max(res$nsnp)
		lam$minp[i] <- min(res$pval)
		lam$nsig[i] <- sum(res$fdr < 0.05)
		lam$lambda[i] <- lambda(res$pval)[1]
		lam$proppos[i] <- sum(sign(res$b)==1) / nrow(res)
		lam$pdir[i] <- binom.test(x=sum(sign(res$b) == 1), n=nrow(res), p=0.5)$p.value
	}
}
lam <- merge(lam, tr, by)

save(lam, file="../results/lambda.rdata")
load("../results/lambda.rdata")

table(lam$pdir < 0.05)
hist(lam$proppos)
plot(lambda ~ I(abs(0.5 - proppos)), subset(lam, !chr6))
summary(lm(lambda ~ I(abs(0.5 - proppos)), subset(lam, !chr6)))

lama <- subset(lam, chr6)
lamb <- subset(lam , !chr6)
plot(lama$proppos, lamb$proppos)
library(ggplot2)
ggplot(lam, aes(x=lambda, y=nsig)) +
geom_point(aes(colour=chr6))

# need to exclude chr6 - leading to massive numbers of hits
summary(lm(lambda ~ I(abs(0.5 - proppos)), lamb))
plot(nsig ~ lambda, lamb)

subset(lamb, lambda > 1.05)
lamb$phen <- do.call(c, lapply(as.character(lamb$trait), function(x) strsplit(x, split="\\|")[[1]][1])) %>% gsub(" $", "", .)
lamb <- lamb %>% arrange(desc(lambda)) %>% subset(!duplicated(phen))
lamb$phen <- as.factor(lamb$phen)
lamb$phen <- factor(lamb$phen, levels=lamb$phen[order(lamb$lambda)])
ggplot(lamb %>% arrange(lambda), aes(x=phen, y=lambda)) +
geom_point() +
geom_hline(yintercept=1, linetype="dotted") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="Phenotype", y="Inflation factor")
ggsave("../images/lambda.pdf", width=15, height=8)

library(ggrepel)
lamb$imbalance <- abs(lamb$proppos - 0.5)
ggplot(lamb, aes(x=lambda, y=imbalance)) +
geom_point() +
geom_smooth(method="lm") +
geom_label_repel(data=subset(lamb, lambda > 1.05 & imbalance > 0.06), aes(label=phen), size=2) +
labs(x="p-value inflation", y="Enrichment of positive or negative causal effects")
ggsave("../images/lambda_vs_imbalance.pdf")





subset(lamb, lambda>1.05)

lambtop <- subset(lamb, lambda > 1.05)
lambtop

l <- list()
for(i in 1:nrow(lambtop))
{
	message(i)
	load(paste0("../results/out/gwas", lambtop$id[i], ".rdata"))
	res <- res[sample(1:nrow(res)), ]
	res$cut <- cut(1:nrow(res), 500)
	out <- res %>% group_by(cut) %>%
	summarise(lambda=lambda(pval)[1])
	out$id <- lambtop$id[i]
	out$phen <- lambtop$phen[i]
	l[[i]] <- out
}
lambdatop <- bind_rows(l)
save(lambdatop, file="../results/lambda_top.rdata")
load("../results/lambda_top.rdata")

ggplot(lambdatop, aes(x=lambda)) +
geom_histogram() +
geom_vline(data=lambtop, aes(xintercept=lambda), colour="red") +
geom_vline(xintercept=1, linetype="dotted") +
facet_wrap(~ phen) +
theme(strip.text=element_text(size=6))
ggsave("../images/lambda_top.pdf")

hist(lam32$lambda, breaks=30)


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

