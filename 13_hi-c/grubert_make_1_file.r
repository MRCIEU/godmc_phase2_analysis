library(data.table)

# there should be 40041300-130(header) lines in concatenated file
files <- list.files(pattern = ".tsv")
temp <- lapply(files, fread)
hic <- rbindlist( temp )
str(hic)

save(hic, file="Grubert_HiC_0.4_clean.Rdata")