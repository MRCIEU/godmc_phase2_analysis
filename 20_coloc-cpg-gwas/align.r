library(data.table)
library(tidyr)
library(dplyr)
library(gwasvcf)
set_bcftools()
stopifnot(check_bcftools())

args <- commandArgs(T)
infile <- args[1]
ref <- args[2]
outfile <- args[3]

# infile <- "16_100.txt.gz"
# ref <- "/mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/data/ref/eur.vcf.gz"

a <- data.table::fread(infile, he=T) %>%
	tidyr::separate(MarkerName, c("chr", "pos", "rem"), ":") %>%
	tidyr::separate(rem, c("type", "cpg"), "_") %>%
	dplyr::mutate(chr = gsub("chr", "", chr)) %>%
	dplyr::filter(nchar(Allele1) == 1)

reference <- query_gwas(ref, chrompos=paste(a$chr, a$pos, sep=":") %>% unique) %>% VariantAnnotation::expand()
names(reference) <- VariantAnnotation::info(reference)$RS  %>% paste0("rs", .)

temp <- gwasvcf::harmonise(
	chr1 = a$chr, 
	pos1 = as.numeric(a$pos), 
	ref1 = a$Allele2, 
	alt1 = a$Allele1, 
	chr2 = reference %>% SummarizedExperiment::rowRanges() %>% GenomicRanges::seqnames() %>% as.character, 
	pos2 = reference %>% SummarizedExperiment::rowRanges() %>% GenomicRanges::ranges() %>% {.@start}, 
	ref2 = reference %>% SummarizedExperiment::rowRanges() %>% {.$REF} %>% as.character, 
	alt2 = reference %>% SummarizedExperiment::rowRanges() %>% {.$ALT} %>% as.character,
	rsid2 = names(reference),
	indel_recode=TRUE, 
	strand_flip=TRUE
)

a$ref <- temp$ref.y
a$alt <- temp$alt.y
a$decision <- temp$decision

index <- a$decision %in% c(1, 3, 5)

a$Effect[index] <- a$Effect[index] * -1
a$EffectARE[index] <- a$EffectARE[index] * -1
a$Direction[index] <- gsub("\\+", "X", a$Direction[index])
a$Direction[index] <- gsub("-", "+", a$Direction[index])
a$Direction[index] <- gsub("X", "-", a$Direction[index])
a$Freq1[index] <- 1 - a$Freq1[index]
a$dbsnpid <- temp$rsid

x <- a$Allele1[index]
a$Allele1[index] <- a$Allele2[index]
a$Allele2[index] <- x

a$snp <- paste0("chr", paste(a$chr, a$pos, a$type, sep=":"))
a$pos <- as.numeric(a$pos)
a$Pvalue <- as.numeric(a$Pvalue)

save(a, file=outfile)



