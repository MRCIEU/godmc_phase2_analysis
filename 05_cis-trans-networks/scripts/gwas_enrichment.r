library(dplyr)
library(data.table)
library(parallel)


a <- file.path("../data/extract_gwas", list.files(pattern="*enr.rdata", "../data/extract_gwas"))
file.exists(a)

o <- list()
for(i in 1:length(a))
{
	load(a[i])
	res$id <- gsub("../data/extract_gwas/", "", a[i]) %>% gsub("_enr.rdata", "", .)
	o[[i]] <- res
}

gwas_enrichment <- bind_rows(o)
names(gwas_enrichment) <- c("lor", "se", "z", "p", "background", "cluster", "ncase", "ncontrol", "id")

min(gwas_enrichment$p)
group_by(gwas_enrichment, background, is.na(cluster)) %>% summarise(m=min(p))

gwas_enrichment <- group_by(gwas_enrichment, background, is.na(cluster)) %>% mutate(fdr=p.adjust(p, "fdr"))
gwas_enrichment %>% summarise(n=n(), nsig=sum(fdr < 0.05))
subset(gwas_enrichment, fdr < 0.05) %>% as.data.frame

save(gwas_enrichment, file="../results/gwas_enrichment.rdata")


