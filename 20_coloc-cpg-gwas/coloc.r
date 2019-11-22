library(gwasglue)
library(gwasvcf)
library(dplyr)
library(magrittr)
library(coloc)
library(parallel)
library(S4Vectors)
set_bcftools()


args <- commandArgs(T)
infile <- args[1]
outfile <- args[2]

load(infile)
load("../results/16/16_clumped.rdata")
ids <- scan("ids.txt", what="character")
gwas <- tibble(
	id = ids,
	fn = paste0("/mnt/storage/private/mrcieu/research/scratch/IGD/data/public/", id, "/", id, ".vcf.gz")
)
gwas <- gwas[file.exists(gwas$fn),]


b <- subset(a, paste(snp, cpg) %in% paste(clumped$snp, clumped$cpg))
chrompos <- parse_chrompos(paste0(b$chr, ":", b$pos), radius=250000)
values(chrompos)[["cpg"]] <- b$cpg

cpg <- a %>% subset(., cpg %in% b$cpg & !is.na(ref)) %>%
	split(., f = .$cpg) %>%
	lapply(., function(x)
	{
		x %$% create_vcf(chrom=chr, pos=as.numeric(pos), ea=alt, nea=ref, snp=snp, name=cpg[1], ea_af=Freq1, effect=Effect, se=StdErr, pval=Pvalue, n=TotalSampleSize)
	})

chrompos <- chrompos[chrompos$cpg %in% names(cpg)]

l <- list()
for(i in 1:nrow(gwas))
{
	message(gwas$id[i])
	out <- query_gwas(gwas$fn[i], chrompos)
	out <- sort(out[!duplicated(out)])
	l[[gwas$id[i]]] <- lapply(1:length(chrompos), function(j)
	{
		message(j)
		g <- try(gwasvcf_to_coloc(cpg[[values(chrompos)$cpg[[j]]]], out, chrompos[j]))
		try(coloc::coloc.abf(g[[1]], g[[2]])$summary)
	})
}
save(l, chrompos, b, file=outfile)
