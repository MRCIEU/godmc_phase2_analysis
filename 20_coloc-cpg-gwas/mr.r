library(gwasvcf)
library(GenomicRanges)
library(TwoSampleMR)
library(dplyr)
library(magrittr)
library(parallel)
set_bcftools()

args <- commandArgs(T)
i <- as.numeric(args[1])

load("../results/16/16_clumped_aligned.rdata")

ids <- scan("ids.txt", what="character")
gwas <- tibble(
	id = ids,
	fn = paste0("/mnt/storage/private/mrcieu/research/scratch/IGD/data/public/", id, "/", id, ".vcf.gz")
)
gwas <- gwas[file.exists(gwas$fn),]

e <- clumped_aligned %$% tibble(
	SNP = paste0(snpchr, ":", snppos),
	effect_allele.exposure = alt,
	other_allele.exposure = ref,
	beta.exposure = Effect,
	se.exposure = StdErr,
	pval.exposure = pval,
	samplesize.exposure = TotalSampleSize,
	mr_keep.exposure = TRUE,
	id.exposure = cpg,
	exposure = cpg
)

chrompos <- parse_chrompos(paste0(clumped_aligned$snpchr, ":", clumped_aligned$snppos))

mcols(chrompos)$cpg <- clumped_aligned$cpg
mcols(chrompos)$snp <- clumped_aligned$snp
mcols(chrompos)$dbsnpid <- clumped_aligned$dbsnpid
chrompos <- sort(chrompos)

het <- list()
mrivw <- list()
wr <- list()
message(gwas$id[i])
out <- query_gwas(gwas$fn[i], chrompos)
out <- sort(out[!duplicated(out)])
o <- vcf_to_granges(out)
mcols(o)$dbsnpid <- names(o)
o <- as_tibble(o)
o <- o %$% dplyr::tibble(
	SNP=paste0(seqnames, ":", start), 
	other_allele.outcome = REF, 
	effect_allele.outcome = ALT, 
	pval.outcome = 10^-LP,
	samplesize.outcome = SS,
	beta.outcome = ES,
	se.outcome = SE,
	mr_keep.outcome = TRUE,
	id.outcome = gwas$id[i],
	outcome = gwas$id[i])

dat <- inner_join(e, o, by="SNP")
table((dat$effect_allele.exposure == dat$effect_allele.outcome & dat$other_allele.exposure == dat$other_allele.outcome))

index <- dat$other_allele.exposure == dat$effect_allele.outcome & dat$effect_allele.exposure == dat$other_allele.outcome
index[is.na(index)] <- FALSE
dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
dat$mr_keep <- TRUE

mcpg <- table(dat$exposure)
mcpg <- names(mcpg)[mcpg > 1]
datm <- subset(dat, exposure %in% mcpg)
mrivw <- mr(datm, method_list=c("mr_ivw"))
het <- mr_heterogeneity(datm, method_list="mr_ivw")
wr <- dat %>%
	mutate(id.exposure = paste(SNP, exposure), exposure = paste(SNP, exposure)) %>% mr(., method_list=c("mr_ivw", "mr_wald_ratio"))

save(mrivw, het, wr, file=paste0("results/scratch/mr", i, ".rdata"))

