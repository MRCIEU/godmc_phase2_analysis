library(tidyverse)
library(TwoSampleMR)

a <- read_tsv("17_92.txt.gz")
load("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")

bmi <- subset(gwas, grepl("Body mass", exposure) & data_source.exposure == "mrbase")

b <- separate(a, MarkerName, c("snp", "cpg"), "_")

temp <- subset(b, snp %in% bmi$id)
temp <- format_data(temp, type="outcome",
	phenotype_col="cpg",
	snp_col="snp",
	beta_col="Effect",
	se_col="StdErr",
	effect_allele_col="Allele1",
	other_allele_col="Allele2",
	eaf_col="Freq1",
	pval_col="P-value"
)

bmi$SNP <- tolower(bmi$id)
bmi$exposure <- bmi$trait
bmi <- subset(bmi, select=-c(id, trait))
bmi$id.exposure <- 1

dat <- harmonise_data(bmi, temp, action=1)
m <- mr(dat, metho=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
table(m$pval < 0.001)
min(m$pval)
min(p.adjust(m$pval))


