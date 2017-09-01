library(dplyr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(parallel)

setwd("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/scripts/godmc_phase2_analysis/06_mr-gwas-cpg/scripts/")

load("/mnt/storage/home/gh13047/repo/godmc_phase2_analysis/06_mr-gwas-cpg/results/tophits.rdata")

load("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")
a <- subset(gwas, mr_keep.exposure == TRUE)
# a <- subset(a, exposure %in% res$exposure)

res$code <- paste(res$exposure, res$outcome)

param <- expand.grid(
	gwas=unique(a$exposure),
	chunk=1:300
)

args <- commandArgs(T)

jid <- as.numeric(args[1])
num <- as.numeric(args[2])

first <- (jid - 1) * num + 1
last <- min(jid * num, nrow(param))

param <- param[first:last, ]

out <- paste0("../results/out_followup", jid, ".rdata")
outtemp <- paste0(out, ".temp")

if(file.exists(out)) q()

if(file.exists(outtemp))
{
	load(outtemp)
	i1 <- length(m) + 1
} else {
	i1 <- 1
}

curr <- "0"
m <- list()
for(i in i1:nrow(param))
{
	message(i, " of ", nrow(param))
	if(param$gwas[i] %in% res$exposure)
	{
		fn <- paste0("zcat ../../results/17/17_", param$chunk[i], ".txt.gz")
		if(fn != curr)
		{
			message("Reading")
			b <- fread(fn)
			b <- separate(b, MarkerName, c("snp", "cpg"), "_")
			b <- subset(b, cpg %in% res$outcome)
			curr <- fn
		} else {
			message("Same cpg file")
		}
		exposure <- subset(a, exposure == param$gwas[i])
		exposure$SNP <- exposure$id
		exposure$exposure <- param$gwas[i]
		exposure <- subset(exposure, select=-c(id, trait))

		outcome <- subset(b, snp %in% exposure$SNP)
		if(nrow(outcome) > 0)
		{
			message("Formatting")
			outcome <- format_data(outcome, type="outcome",
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
			exposure$SNP <- tolower(exposure$SNP)
			message("Harmonising")
			dat <- suppressMessages(harmonise_data(exposure, outcome, action=2))
			dat$code <- paste(dat$exposure, dat$outcome)
			dat <- subset(dat, code %in% res$code)
			if(nrow(dat) > 0)
			{
				message("Analysing")
				m[[i]] <- suppressMessages(mr(dat))
				message("Saving")
				save(m, dat, file=outtemp)
			}
		}
	}
}

res <- bind_rows(m)

save(res, file=out)
unlink(outtemp)

