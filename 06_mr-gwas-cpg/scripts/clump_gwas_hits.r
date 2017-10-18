library(TwoSampleMR)
library(dplyr)
library(biomaRt)
library(GenomicRanges)

load("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")
a <- subset(gwas, mr_keep.exposure == TRUE)

a1 <- subset(a, data_source.exposure == "mrbase")
a2 <- subset(a, data_source.exposure != "mrbase")

sp <- split(unique(a2$exposure), 1:100)

l <- list()
for(i in 1:100)
{
	l[[i]] <- clump_data(subset(a, exposure %in% sp[[i]]))
}

a2 <- bind_rows(l) %>% filter(!duplicated(paste(SNP, id.exposure)))
a <- rbind(a1, a2)

load("../../data/misc/ensembl_gene_list.rdata")
load("../../05_cis-trans-networks/data/snpcontrolsets_selection.rdata")
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies))
names(ldinfo) <- c("snp", "snplow", "snphigh", "snpproxies")

temp <- subset(a, !duplicated(id))
temp <- inner_join(temp, ldinfo, c("id"="snp"))
temp <- tidyr::separate(temp, id, c("snpchr", "snppos", "snptype"), sep=":", rem=FALSE)
snps <- GRanges(seqnames=as.character(temp$snpchr), ranges=IRanges(temp$snplow, temp$snphigh), strand="*")
names(snps) <- temp$id

b <- findOverlaps(snps, genes)
snp_gene <- data_frame(snp=names(snps)[b@queryHits], snp_hgnc_symbol=mcols(genes)[b@subjectHits,2]) %>%
	filter(!duplicated(snp))
a <- merge(a, snp_gene, by.x="id", by.y="snp", all.x=TRUE)
a$snp_hgnc_symbol[a$snp_hgnc_symbol == ""] <- NA
index <- is.na(a$gene.exposure)
a$gene.exposure[index] <- a$snp_hgnc_symbol[index]
index <- is.na(a$gene.exposure)

save(a, file="../data/snps_gwas.rdata")
