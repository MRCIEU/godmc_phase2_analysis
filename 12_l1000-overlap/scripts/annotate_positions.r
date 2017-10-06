library(biomaRt)
library(dplyr)
library(GenomicRanges)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
genes <- lapply(1:23, function(x)
{
	getBM(
		c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position","end_position"),
		filters = c("chromosome_name"), 
		values = list(x),
		mart = ensembl)
})

genes <- bind_rows(genes)
save(genes, file="../../data/misc/ensembl_gene_list.rdata")

load("../../results/16/16_clumped.rdata")

load("../../05_cis-trans-networks/data/snpcontrolsets_selection.rdata")
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))


temp <- subset(clumped, !duplicated(snp))
snpgr <- 
