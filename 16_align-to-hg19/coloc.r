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
cores <- as.numeric(args[3])

# load("../results/16/16_801_aligned.rdata")
load(infile)
load("../results/16/16_clumped.rdata")

gwas <- tibble(
	id = c("GCST004599", "GCST004600", "GCST004601", "GCST004602", "GCST004603", "GCST004604", "GCST004605", "GCST004606", "GCST004607", "GCST004608", "GCST004609", "GCST004610", "GCST004611", "GCST004612", "GCST004613", "GCST004614", "GCST004615", "GCST004616", "GCST004617", "GCST004618", "GCST004619", "GCST004620", "GCST004621", "GCST004622", "GCST004623", "GCST004624", "GCST004625", "GCST004626", "GCST004627", "GCST004628", "GCST004629", "GCST004630", "GCST004631", "GCST004632", "GCST004633", "GCST004634"),
	fn = paste0("/mnt/storage/private/mrcieu/research/scratch/IGD/data/public/EBI-a-", id, "/EBI-a-", id, ".vcf.gz")
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
	l[[gwas$id[i]]] <- mclapply(1:length(chrompos), function(j)
	{
		message(j)
		g <- vcflist_overlaps(list(cpg[[values(chrompos)$cpg[[j]]]], out), chrompos[j])
		try(perform_coloc(g[[1]], g[[2]])$summary)
	}, mc.cores=cores)
}
save(l, chrompos, b, file=outfile)

