
# BMI, HDL, LDL, trigs
# 

library(dplyr)

ewas <- read.table("../data/EWAS_Catalog_20-02-2018.txt.gz", he=T, sep="\t", stringsAsFactors=FALSE)
ewas$code <- paste(ewas$PMID, ewas$Trait)
b <- group_by(subset(ewas, P < 1e-7), code) %>% summarise(n=n()) %>% arrange(desc(n))

traits2 <- read.csv("../../data/gwas/00info.csv")
traits2$id_10 <- 1:nrow(traits2)
load("../data/snps_gwas.rdata")

library(TwoSampleMR)
ao <- available_outcomes()
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
data.frame(substr(as.character(a1$exposure), start=1, stop=50), a1$n)
remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)
a1 <- a0[-remove_ids,]



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



