library(TwoSampleMR)
ao <- available_outcomes()
library(magrittr)
library(dplyr)

load("~/repo/godmc-database/neo4j/data/trait_id_master.rdata")
load("../results/gwas_enrichment.rdata")
b <- subset(ao, access != "developer" & id %in% subset(master, !is.na(id_06))$id_mrb)
uid <- unique(gwas_enrichment$id)
uid <- uid[!uid %in% b$id]


b0 <- b %$% 
	data_frame(trait, author, pmid, sample_size, ncase, ncontrol, subcategory, mr_outcome=TRUE)
b1 <- subset(ao, id %in% uid) %$%
	data_frame(trait, author, pmid, sample_size, ncase, ncontrol, subcategory, mr_outcome=FALSE)

b0$mr_exposure <- b0$mr_outcome

b01 <- bind_rows(b0, b1)

remove <- c(
	grep("Father", b01$trait),
	grep("Mother", b01$trait),
	grep("Parent", b01$trait),
	grep("survival", b01$trait),
	grep("Difference", b01$trait),
	grep("Subjective", b01$trait),
	grep("Melanoma", b01$trait)
)
remove
b01$mr_exposure[remove] <- FALSE
dim(b01)
write.csv(b01, "../results/trait_list.csv")


