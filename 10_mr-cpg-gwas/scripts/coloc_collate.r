library(data.table)
library(dplyr)

fn <- paste0("../results/coloc/coloc", 1:158, ".rdata")
l <- list()
for(i in fn)
{
	if(file.exists(i))
	{
		load(i)
		l[[i]] <- res
	} else {
		message(i, " missing")
	}
}

res <- bind_rows(l)
res <- subset(res, !is.na(H4))
ext <- fread("zcat ../results/filtered_gwas_mqtl.txt")
names(ext) <- c("outcome", "snp", "ea", "oa", "eaf", "b", "se", "p", "n")
names(res)[1] <- "snp"
ext$outcome <- as.character(ext$outcome)

res <- inner_join(res, ext, by=c("snp", "outcome"))

sum(res$p < 1e-10)
sum(res$H4 > 0.8)
sum(res$H4 > 0.8 & res$p < 1e-10)
sig <- res[res$H4 > 0.8 & res$p < 1e-10,]
save(res, file="../results/cpg_trait_coloc.rdata")



library(TwoSampleMR)
ao <- available_outcomes()

ao <- subset(ao, id %in% res$outcome)
sig <- res[res$H4 > 0.8 & res$p < 1e-10,]
sig <- subset(sig, !duplicated(paste(exposure, outcome)))

sig <- merge(sig, subset(ao, select=c(id, trait, category, subcategory)), by.x="outcome", by.y="id")

sigm <- subset(sig, category != "Metabolites")
temp <- table(sigm$exposure) %>% sort(decreasing=TRUE)

subset(sigm, exposure == names(temp)[1])$trait
subset(sigm, exposure == names(temp)[2])$trait
subset(sigm, exposure == "cg26489994")$trait
subset(sigm, exposure == names(temp)[3])
subset(sigm, exposure == names(temp)[4])

