library(dplyr)
library(data.table)
library(TwoSampleMR)

load("../results/to_extractt.rdata")
load("../results/to_extractc.rdata")
load("../data/snp_1kg.rdata")

to_extractt <- inner_join(to_extractt, snp_1kg, by=c("transsnp"="snp"))
id <- unique(to_extractt$outcome)
l1 <- list()
for(i in 1:length(id))
{
	message(i, " of ", length(id))
	temp <- subset(to_extractt, outcome == id[i])
	l1[[i]] <- extract_outcome_data(temp$V2, id[i])
}

trans_extract <- bind_rows(l1)
temp <- subset(to_extractt, !duplicated(transsnp), select=c(transsnp, V2))
names(temp) <- c("transsnp", "SNP")
trans_extract <- merge(trans_extract, temp, by="SNP")


to_extractc <- inner_join(to_extractc, snp_1kg, by=c("transsnp"="snp"))
id <- unique(to_extractc$outcome)
l2 <- list()
for(i in 1:length(id))
{
	message(i, " of ", length(id))
	temp <- subset(to_extractc, outcome == id[i])
	l2[[i]] <- extract_outcome_data(temp$V2, id[i])
}

cis_extract <- bind_rows(l2)
temp <- subset(to_extractc, !duplicated(transsnp), select=c(transsnp, V2))
names(temp) <- c("transsnp", "SNP")
cis_extract <- merge(cis_extract, temp, by="SNP")


save(trans_extract, file="../results/trans_extract.rdata")
save(cis_extract, file="../results/cis_extract.rdata")
