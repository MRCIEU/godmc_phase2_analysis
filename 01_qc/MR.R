#~/Documents/MRBase
library(TwoSampleMR)
toggle_dev("elastic")

snps<- scan("~/Documents/MRBase/rsids_Mstat.txt",what="character")
g<-grep("duplicate",snps)
snps<-snps[-g]
#snps<- scan("rsids_Mstat.txt",what="character")
snps
ids<-read.table("ids.txt") #254
ao <- available_outcomes()
ids <- ids[which(ids$V1 %in% ao$id),] #230

b <- extract_outcome_data(snps, 2)
b <- extract_outcome_data(snps, ids)
#b <- extract_outcome_data(snps, ids$V1,proxies=FALSE)
head(b)
mean(b$pval.outcome)
hist(b$pval.outcome)
min(b$pval.outcome)
hist(b$pval.outcome, break=100)
hist(b$pval.outcome, breaks=100)
b[which(b$pval.outcome<5e-8),]