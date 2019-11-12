#library(TwoSampleMR)
#ao <- available_outcomes()
library(magrittr)
library(dplyr)

#save(ao, file="../data/outcomes.rdata")


load("../data/outcomes.rdata")
load("~/repo/godmc-database/neo4j/data/trait_id_master.rdata")
load("../results/gwas_enrichment.rdata")
b <- subset(ao, access != "developer" & id %in% subset(master, !is.na(id_06))$id_mrb)
uid <- unique(gwas_enrichment$id)
uid <- uid[!uid %in% b$id]


b0 <- b %$% 
	data_frame(trait, author, pmid, sample_size, ncase, ncontrol, subcategory, id, mr_outcome=TRUE)
b1 <- subset(ao, id %in% uid) %$%
	data_frame(trait, author, pmid, sample_size, ncase, ncontrol, subcategory, id, mr_outcome=FALSE)

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

g<-grep("Type 1 diabetes",b01$trait)
b01[g,"author"]<-"Barrett JC/ Cooper JD"
b01[g,"pmid"]<-"19430480/18978792"

#b01<-b01[which(!is.na(b01$mr_exposure)&!is.na(b01$mr_outcome)),]
#b01<-b01[which(b01$mr_outcome=="TRUE"| b01$mr_exposure=="TRUE"),]

community<-read.csv("../results/gwas_enrichment.csv")
w<-which(b01$id%in%community$mrbase_id)
b01$community<-"FALSE"
b01$community[w]<-"TRUE"
length(unique(community$mrbase_id)) #133
b01<-b01[which(b01$mr_outcome=="TRUE"| b01$mr_exposure=="TRUE"|b01$community=="TRUE"),]
length(which(b01$mr_outcome=="TRUE"| b01$mr_exposure=="TRUE"|b01$community=="TRUE")) #144
write.csv(b01, "../results/trait_list.csv")

