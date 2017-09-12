# 1. Get surrounding SNP-CpG associations for putative SNP
# 2. Get those SNPs for SNP-trait associations
# 3. Get correlation matrix for remaining SNPs
# 4. Run coloc

library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(coloc)
options(dplyr.show_progress = FALSE)


convert_to_rs <- function(snp, snp_1kg=snp_1kg)
{
	index <- match(snp, snp_1kg$snp)
	return(snp_1kg$V2[index])
}

convert_to_chrpos <- function(snp, snp_1kg=snp_1kg)
{
	index <- match(snp, snp_1kg$V2)
	return(snp_1kg$snp[index])
}

run_coloc <- function(dat, type1, type2)
{

	# H0 (no causal variant), 
	# H1 (causal variant for trait 1 only), 
	# H2 (causal variant for trait 2 only), 
	# H3 (two distinct causal variants)
	# H4 (one common causal variant)

	require(coloc)
	dataset1 <- list(
		pvalues=dat$pval.exposure,
		N=dat$samplesize.exposure,
		beta=dat$beta.exposure,
		varbeta=dat$se.exposure^2,
		snp=dat$SNP,
		type=type1,
		MAF=dat$eaf.exposure
	)
	dataset2 <- list(
		pvalues=dat$pval.outcome,
		N=dat$samplesize.outcome,
		beta=dat$beta.outcome,
		varbeta=dat$se.outcome^2,
		snp=dat$SNP,
		type=type2,
		MAF=dat$eaf.outcome
	)
	dataset2$MAF[is.na(dataset2$MAF)] <- dataset1$MAF[is.na(dataset2$MAF)]
	a <- suppressMessages(coloc.abf(dataset1, dataset2))
	return(data.frame(
		exposure=dat$exposure[1],
		outcome=dat$outcome[1],
		nsnp=a$summary[1],
		H0=a$summary[2],
		H1=a$summary[3],
		H2=a$summary[4],
		H3=a$summary[5],
		H4=a$summary[6]
	))
}

## Setup data


## Parallel

args <- commandArgs(T)
jid <- as.numeric(args[1])
num <- as.numeric(args[2])

out <- paste0("../results/coloc", jid, ".rdata")

if(file.exists(out)) q()

load("../data/coloc_data.rdata")
snplist <- unique(filtered_gwas_mqtl$V2)

first <- (jid - 1) * num + 1
last <- min(jid * num, length(snplist))
snplist <- snplist[first:last]

i <- 0
res <- filter(filtered_gwas_mqtl, V2 %in% snplist) %>%
	group_by(V2) %>%
	do(
	{
		i <- i + 1
		selected_snp <- .
		candidate_cpgs <- subset(clumped, snp2 == selected_snp$V2[1])
		ncpg <- nrow(candidate_cpgs)
		ntrait <- nrow(selected_snp)
		for(cpgnum in 1:ncpg)
		{
			unclumped <- fread(paste0("zcat ../../results/16/16_", candidate_cpgs$chunk[cpgnum], ".txt.gz")) %>%
				separate(MarkerName, c("snp", "cpg"), "_") %>%
				filter(cpg == candidate_cpgs$cpg[cpgnum], snp %in% snp_1kg$snp)
			unclumped$snp2 <- convert_to_rs(unclumped$snp, snp_1kg)
			unclumped <- suppressMessages(format_data(
				unclumped,
				type="exposure",
				snp_col="snp2",
				phenotype_col="cpg",
				effect_allele_col="Allele1",
				other_allele_col="Allele2",
				eaf_col="Freq1",
				beta_col="Effect",
				se_col="StdErr",
				pval_col="Pvalue",
				samplesize_col="TotalSampleSize"
			))
			l <- list()
			j <- 1
			for(traitnum in 1:ntrait)
			{
				fn <- paste0("../scratch/filtered_gwas_", selected_snp$V1[traitnum], ".txt")
				trait <- fread(fn) %>%
					filter(snp %in% unclumped$SNP)
				message(i, " : ", ncpg, " : ", ntrait, " : ", selected_snp$V2[1], " : ", candidate_cpgs$cpg[cpgnum], " : ", selected_snp$V1[traitnum])
				if(nrow(trait) > 0)
				{				
					trait <- suppressMessages(format_data(
						trait,
						type="outcome",
						snp_col="snp",
						eaf_col="effect_allele_freq",
						pval_col="p",
						samplesize_col="n"
					))
					trait$outcome <- selected_snp$V1[traitnum]
					dat <- suppressMessages(harmonise_data(unclumped, trait))
					res <- try(run_coloc(dat, "quant", ifelse(selected_snp$cc[traitnum], "cc", "quant")))
					if(class(res) == 'try-error')
					{
						l[[j]] <- NULL
					} else {
						l[[j]] <- res
					}
				}
				j <- j+1
			}
		}
		return(bind_rows(l))
	})

save(res, file=out)
