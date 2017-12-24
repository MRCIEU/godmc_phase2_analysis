suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))

remove_bad_snps <- function(x, bim)
{
	indel <- subset(x, grepl("INDEL", snp))
	snps <- subset(x, grepl("SNP", snp))
	bim <- subset(bim, snp %in% snps$snp)
	snps <- subset(snps, snp %in% bim$snp)
	index <- match(snps$snp, bim$snp)
	bim <- bim[index, ]
	stopifnot(all(bim$snp == snps$snp))	

	palin <- (
		(snps$Allele1 == "a" & snps$Allele2 == "t") |
		(snps$Allele1 == "t" & snps$Allele2 == "a") |
		(snps$Allele1 == "g" & snps$Allele2 == "c") |
		(snps$Allele1 == "g" & snps$Allele2 == "c")) &
		(snps$Freq1 > 0.45 & snps$Freq1 < 0.55)


	flipped <- (snps$Allele1 != bim$a1 & snps$Allele2 != bim$a2) &
		(snps$Allele2 != bim$a1 & snps$Allele1 != bim$a2)

	message(nrow(indel), " indels")
	message(nrow(snps), " snps")
	message(sum(palin), " palindromes to exclude")
	message(sum(flipped), " flipped to exclude")
	snps <- snps[!(palin | flipped), ]
	return(rbind(indel, snps))
}

cleaning_message <- function(...)
{
	message("CLEANING MESSAGE: ", paste(..., collapse=""))
}

i <- as.numeric(commandArgs(T)[1])
out <- paste0("../results/16/16_cleaned_", i, ".rdata")
# if(file.exists(out)) q()

# Get cpg positions
load("../data/misc/cpg_pos.rdata")
# load("../data/ref/alleles.rdata")
# bim$a1 <- tolower(bim$a1)
# bim$a2 <- tolower(bim$a2)

# Only retain zhou CpGs
zhou <- scan("../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())
twinsuk <- fread("../../godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt")$ILMID
badsnps <- scan("../04_conditional_16/badsnps.txt", what=character())

cis_radius <- 1000000

message(i)
a <- fread(paste0("zcat ../results/16/16_", i, ".txt.gz"))
a$Pvalue <- as.numeric(a$Pvalue)
a$HetPVal <- as.numeric(a$HetPVal)
a$PvalueARE <- as.numeric(a$PvalueARE)
a$PvalueMRE <- as.numeric(a$PvalueMRE)

a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
a$snp2 <- a$snp
a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
a$snppos <- as.numeric(a$snppos)
a <- inner_join(a, cpgpos, by=c("cpg"))
a$cis <- FALSE
a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE
# a <- remove_bad_snps(a, bim)

cleaning_message(nrow(a), " total associations")
cleaning_message(length(unique(a$snp)), " unique SNPs")
cleaning_message(length(unique(a$cpg)), " unique CPGs")

# remove indels of same length
ind <- a$Allele1 %in% c("a", "c", "t", "g", "d", "i") & a$Allele2 %in% c("a", "c", "t", "g", "d", "i")
a <- a[ind, ]
cleaning_message(sum(!ind), " INDELs of same length being removed")

ind <- a$snp %in% badsnps
cleaning_message(sum(ind), " SNPs being removed for being multi-allelic")
a <- a[!ind, ]

# Remove results with fewer than 5 studies or fewer than 5k samples
nstudies <- nchar(a$Direction[1])
nst <- nstudies - str_count(a$Direction, "\\?")
rem <- nst < 5 | a$TotalSampleSize < 5000
cleaning_message(sum(rem), " associations with fewer than 5 studies or 5000 samples removed")

res <- a[!rem, ]

ind <- res$cpg %in% zhou
cleaning_message(sum(!ind), " associations removed for being bad CpGs (Zhou)")

res <- subset(res, ind)

ind <- ! res$cpg %in% twinsuk
cleaning_message(sum(!ind), " further associations removed for being bad CpGs (TwinsUK)")

res <- subset(res, ind)

cleaning_message(nrow(res), " associations remaining")
save(res, file=out)
