library(TwoSampleMR)
ao <- available_outcomes()
library(magrittr)
library(dplyr)

load("../data/outcomes.RData")
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

b01$subcategory[b01$subcategory=="Hemodynamic"] <- "Haematological"
b01$subcategory[b01$subcategory=="Haemotological"] <- "Haematological"
b01$subcategory[b01$subcategory=="Immune system"] <- "Autoimmune / inflammatory"
b01$subcategory[b01$subcategory=="Diabetes"] <- "Glycemic"
b01$subcategory[b01$subcategory=="Biomarker"] <- "Other"
b01$subcategory[b01$subcategory=="Protein"] <- "Other"
b01$subcategory[b01$subcategory=="Hormone"] <- "Other"
b01$subcategory[b01$subcategory=="Reproductive aging"] <- "Aging"
b01$subcategory[b01$subcategory=="Lung disease"] <- "Other"
b01$subcategory[b01$subcategory=="Autoimmune / inflammatory"] <- "Immune"
b01$subcategory[b01$subcategory=="Psychiatric / neurological"] <- "Neurological"
b01$subcategory[is.na(b01$subcategory)] <- "Kidney"

write.csv(b01, "../results/trait_list.csv")


##

head(gwas_enrichment)
gwas_enrichment <- inner_join(gwas_enrichment, subset(b, select=c(id, trait, pmid, subcategory))) %>% ungroup() %>%
dplyr::select(community=cluster,trait,pmid, subcategory,  background, ncase, ncontrol, lor, se, z, pval=p, fdr) %>% arrange(pval)
gwas_enrichment$community <- as.character(gwas_enrichment$community)
gwas_enrichment$community[is.na(gwas_enrichment$community)] <- "All"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Hemodynamic"] <- "Haematological"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Haemotological"] <- "Haematological"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Immune system"] <- "Autoimmune / inflammatory"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Diabetes"] <- "Glycemic"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Biomarker"] <- "Other"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Protein"] <- "Other"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Hormone"] <- "Other"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Reproductive aging"] <- "Aging"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Lung disease"] <- "Other"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Autoimmune / inflammatory"] <- "Immune"
gwas_enrichment$subcategory[gwas_enrichment$subcategory=="Psychiatric / neurological"] <- "Neurological"
gwas_enrichment$subcategory[is.na(gwas_enrichment$subcategory)] <- "Kidney"

write.csv(gwas_enrichment, file="../results/gwas_enrichment.csv")
