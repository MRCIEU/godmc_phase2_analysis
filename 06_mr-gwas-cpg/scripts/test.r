library(tidyverse)
library(TwoSampleMR)
library(parallel)

load("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")
bmi <- subset(gwas, grepl("Body mass", exposure) & data_source.exposure == "mrbase")
bmisnp <- bmi$id
bmi$SNP <- tolower(bmi$id)
bmi$exposure <- bmi$trait
bmi <- subset(bmi, select=-c(id, trait))
bmi$id.exposure <- 1

m <- list()
for(i in 1:300)
{
	message(i)
	b <- try(fread(paste0("zcat ../../results/17/17_", i, ".txt.gz")))
	if(nrow(b) > 0 )
	{
		b <- separate(b, MarkerName, c("snp", "cpg"), "_")

		temp <- subset(b, snp %in% bmisnp)
		if(nrow(temp) > 0)
		{
			temp <- format_data(temp, type="outcome",
				phenotype_col="cpg",
				snp_col="snp",
				beta_col="Effect",
				se_col="StdErr",
				effect_allele_col="Allele1",
				other_allele_col="Allele2",
				eaf_col="Freq1",
				pval_col="P-value",
				samplesize_col="TotalSampleSize"
			)



			dat <- harmonise_data(bmi, temp, action=1)
			m[[i]] <- mr(dat, metho=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

		}
	}
}


