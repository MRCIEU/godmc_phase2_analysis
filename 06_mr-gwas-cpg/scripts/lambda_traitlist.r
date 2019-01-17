library(TwoSampleMR)
ao <- available_outcomes()
# ao <- subset(ao, access != "developer")
library(dplyr)

ewas <- read.table("../data/EWAS_Catalog_20-02-2018.txt.gz", he=T, sep="\t", stringsAsFactors=FALSE)
ewas$code <- paste(ewas$PMID, ewas$Trait)
b <- group_by(subset(ewas, P < 1e-7), code) %>% summarise(n=n()) %>% arrange(desc(n))

traits2 <- read.csv("../../data/gwas/00info.csv")
traits2$id_10 <- 1:nrow(traits2)
load("../data/snps_gwas.rdata")

ao$name <- paste(ao$trait, "||", ao$consortium, "||", ao$year, "||", ao$unit)
traits <- unique(subset(a, data_source.exposure == "mrbase")$exposure)
# traits <- traits[traits %in% ao$name]
table(traits %in% ao$name)
rename <- traits[which(!traits %in% ao$name)]

a$exposure_orig <- a$exposure
a$exposure[a$exposure == rename[1]] <- "Serum creatinine (eGFRcrea) || CKDGen || 2015 || log ml/min/1.73 m^2"
a$exposure[a$exposure == rename[2]] <- "Serum cystatin C (eGFRcys) || CKDGen || 2015 || log ml/min/1.73 m^2"
a$exposure[a$exposure == rename[3]] <- "Father's age at death || UK Biobank || 2016 || SD"
a$exposure[a$exposure == rename[4]] <- "Mother's age at death || UK Biobank || 2016 || SD"

a <- subset(a, !exposure %in% subset(ao, access == "developer")$name)


a0 <- filter(a, data_source.exposure == "mrbase") %>% group_by(id.exposure, exposure) %>% summarise(n=n()) %>% arrange(exposure, n)
a0$id_06 <- 1:nrow(a0)

remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)

a0 <- subset(a0, !id_06 %in% remove_ids)

table(a0$exposure %in% ao$name)

temp <- merge(subset(traits2, select=c(id, id_10)), ao, by="id")
table(a0$exposure %in% temp$name)

a0 <- merge(subset(a0, select=c(exposure, id_06)), temp, by.x="exposure", by.y="name", all=TRUE)

temp <- subset(a0, !is.na(id_06) | !is.na(id_10))

# Some traits are duplicated, look by eye to remove the out of date ones
remove_ids <- c(148,147,141,142,137,138,136,134,124,114,111,106,85,86,87,63,64,56,57,33,34,30,20,21,22,18,14,7)
a1 <- a0[-remove_ids,]
data.frame(substr(as.character(a1$exposure), start=1, stop=50), a1$n)


# To remove from trait -> SNP for being unreliable:

# Father's age at death -Pilling LC-PUBMED ID:27015805
# Mother's age at death -Pilling LC-PUBMED ID:27015805
# Parents' age at death -Pilling LC-PUBMED ID:27015805
# Top 1 % survival -Pilling LC-PUBMED ID:27015805
# Difference in height between adolescence and adulthood-Cousminer DL -PUBMED ID:23449627
# Difference in height between childhood and adulthood-Cousminer DL -PUBMED ID:23449627
# Well being-Okbay A-PUBMED ID 27089181
# Melanoma-Amos CI-PUBMED ID 21926416



remove <- c(
	grep("Father", traits),
	grep("Mother", traits),
	grep("Parent", traits),
	grep("survival", traits),
	grep("Difference", traits),
	grep("Subjective", traits),
	grep("Melanoma", traits)
)
remove
traits <- traits[-remove]

dev <- subset(ao, access=="developer")$name
traits <- traits[!traits %in% dev]
save(a0, a1, traits, file="../data/traitlist.rdata")
