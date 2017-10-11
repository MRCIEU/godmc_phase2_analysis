library(TwoSampleMR)
library(dplyr)
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
save(a, file="../data/snps_gwas.rdata")

