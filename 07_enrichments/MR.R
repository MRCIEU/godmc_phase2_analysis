library(TwoSampleMR)
toggle_dev("elastic")
a <- available_outcomes()
snps<- scan("~/Documents/MRBase/rsid.SDS.cisamb.txt",what="character")
snps

ids<-read.table("ids.txt") #254
ao <- available_outcomes()
ids <- ids[which(ids$V1 %in% ao$id),] #230

#l <- list()
#for(i in 1:length(ids))
#{
#    l[[i]] <- extract_outcome_data(snps, ids)
#}

#l <- do.call(rbind, l)
#b <- extract_outcome_data(snps, ids)

b <- extract_outcome_data(snps, ids,proxies=FALSE)
head(b)
mean(b$pval.outcome)
hist(b$pval.outcome)
min(b$pval.outcome)
hist(b$pval.outcome, break=100)
hist(b$pval.outcome, breaks=100)
b2<-b[which(b$pval.outcome<5e-8),]

write.table(b2,"sds.cis.amb.results.txt",sep="\t",col.names=T,row.names=F,quote=F)
