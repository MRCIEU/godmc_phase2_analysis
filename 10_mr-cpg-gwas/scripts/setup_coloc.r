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

filtered_gwas_mqtl <- fread("zcat ../results/filtered_gwas_mqtl.txt.gz")
ao <- available_outcomes()
ao <- subset(ao, id %in% filtered_gwas_mqtl$V1, select=c(id, unit))
ao$cc <- ao$unit == "log odds"
filtered_gwas_mqtl$V1 <- as.character(filtered_gwas_mqtl$V1)
filtered_gwas_mqtl <- merge(filtered_gwas_mqtl, subset(ao, select=c(id, cc)), by.x="V1", by.y="id")
load("../../results/16/16_clumped.rdata")
snp_1kg <- fread("../data/eur.bim.orig")
snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)
i1 <- snp_1kg$c1 == 1 & snp_1kg$c2 == 1
snp_1kg$snp[i1] <- paste0("chr", snp_1kg$V1[i1], ":", snp_1kg$V4[i1], ":SNP")
snp_1kg$snp[!i1] <- paste0("chr", snp_1kg$V1[!i1], ":", snp_1kg$V4[!i1], ":INDEL")
index <- match(clumped$snp, snp_1kg$snp)
stopifnot(all(clumped$snp == snp_1kg$snp[index]))
clumped$snp2 <- snp_1kg$V2[index]

save(snp_1kg, clumped, filtered_gwas_mqtl, file="../data/coloc_data.rdata")

