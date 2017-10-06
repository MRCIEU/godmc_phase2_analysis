library(biomaRt)
library(dplyr)
library(GenomicRanges)

load("../../results/16/16_clumped.rdata")

# Generate GRanges object of SNPs
# Use min/max from LD proxies
load("../../05_cis-trans-networks/data/snpcontrolsets_selection.rdata")
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
names(ldinfo) <- c("snp", "snplow", "snphigh", "snpproxies")

temp <- subset(clumped, !duplicated(snp))
temp <- inner_join(temp, ldinfo, "snp")
snps <- GRanges(seqnames=temp$snpchr, ranges=IRanges(temp$snplow, temp$snphigh), strand="*")
names(snps) <- temp$snp


# Generate GRanges object of CpGs
temp <- subset(clumped, !duplicated(cpg))
cpgs <- GRanges(seqnames=temp$cpgchr, ranges=IRanges(temp$cpgpos, temp$cpgpos), strand="*")
names(cpgs) <- temp$cpg


# Generate GRanges object of genes
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
genes <- GRanges(seqnames=paste0("chr", genes$chromosome_name), ranges=IRanges(genes$start_position, genes$end_position), ensembl_gene_id=genes$ensembl_gene_id, hgnc_symbol=genes$hgnc_symbol)
save(genes, file="../../data/misc/ensembl_gene_list.rdata")


# Find overlaps of SNPs to genes

a <- findOverlaps(snps, genes)
snp_gene <- data_frame(snp=names(snps)[a@queryHits], snp_ensembl_gene_id=mcols(genes)[a@subjectHits,1], snp_hgnc_symbol=mcols(genes)[a@subjectHits,2])
a <- findOverlaps(cpgs, genes)
cpg_gene <- data_frame(cpg=names(cpgs)[a@queryHits], cpg_ensembl_gene_id=mcols(genes)[a@subjectHits,1], cpg_hgnc_symbol=mcols(genes)[a@subjectHits,2])

mqtl <- subset(clumped, select=c(snp, cpg))
mqtl <- merge(mqtl, snp_gene, all=TRUE)
mqtl <- merge(mqtl, cpg_gene, all=TRUE)
mqtl <- as_data_frame(mqtl)
a <- apply(mqtl, 1, function(x) any(is.na(x)))
mqtl <- subset(mqtl, !a)
mqtl <- group_by(mqtl, cpg, snp) %>%
	do({
		x <- .
		y <- expand.grid(
			snp_ensembl_gene_id = unique(x$snp_ensembl_gene_id),
			cpg_ensembl_gene_id = unique(x$cpg_ensembl_gene_id)
		)
		return(y)
	})

save(mqtl, file="../data/mqtl.rdata")
