library(dplyr)

a <- file.path("../data/extract_gwas", list.files(pattern="*rdata", "../data/extract_gwas"))
file.exists(a)

l <- list()
for(i in 1:length(a))
{
	message(i)
	load(a[i])
	l[[i]] <- extracted
}

extracted <- bind_rows(l)
dim(extracted)
save(extracted, file="../data/extract_gwas.rdata")

