library(dplyr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(parallel)

# setwd("/mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/scripts/godmc_phase2_analysis/06_mr-gwas-cpg/scripts/")

load("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")
a <- subset(gwas, mr_keep.exposure == TRUE)


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

out <- paste0("../results/out/out", jid, ".rdata")
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
	fn <- paste0("../../results/17/17_", param$chunk[i], ".txt.gz")
	if(fn != curr)
	{

		message("Extracting")
		snplist <- unique(a$id)
		outlist <- paste0("../scratch/temp", param$chunk[i], "_", jid, ".snplist")
		fnn <- paste0("../scratch/temp", param$chunk[i], "_", jid)
		write.table(snplist, file=outlist, row=F, col=F, qu=F)
		cmd <- paste0("zfgrep -f ", outlist, " ", fn, " > ", fnn)
		system(cmd)
		message("Reading")
		b <- fread(fnn)
		unlink(fnn)
		unlink(outlist)
		nom <- as.character(unlist(read.table(fn, nrows=1)))
		names(b) <- nom
		b$Pvalue <- as.numeric(b$Pvalue)
		b <- separate(b, MarkerName, c("snp", "cpg"), "_")
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
		message("Analysing")
		m[[i]] <- suppressMessages(mr(dat, metho=c("mr_ivw")))
		message("Saving")
		save(m, dat, file=outtemp)
	}
}

res <- bind_rows(m)

save(res, file=out)
unlink(outtemp)
