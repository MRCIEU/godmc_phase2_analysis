library(dplyr)

a <- file.path("../data/extract_gwas", list.files(pattern="*rdata", "../data/extract_gwas"))
file.exists(a)

l <- list()
for(i in 1:length(a))
{
	message(i)
	load(a[i])
	temp <- subset(extracted, duplicated(ID))$ID
	extracted1 <- subset(extracted, !ID %in% temp)
	extracted2 <- subset(extracted, ID %in% temp) %>% filter(is.na(proxy.chrom))
	l[[i]] <- rbind(extracted1, extracted2)
}

extracted <- bind_rows(l)
dim(extracted)
save(extracted, file="../data/extract_gwas.rdata")
