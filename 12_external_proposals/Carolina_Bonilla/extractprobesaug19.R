library(data.table)
library(dplyr)
###

#chunks<-c(1:962)
arguments<-commandArgs(T)
chunk<-as.numeric(arguments[1])

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped<-as.data.frame(clumped)
cpgs<-read.table("cpgs.txt",he=F,stringsAsFactors=F)
chunks<-unique(clumped[which(clumped$cpg%in%cpgs[,1]),"chunk"])

clumped$code <- paste(clumped$snp, clumped$cpg)
#res.all<-data.frame()
r<-data.frame()
for (chunk in 1:length(chunks)){
cat(chunk,"\n")
load(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",chunks[chunk],".rdata",sep=""))
res2<-res[which(res$cpg%in%cpgs[,1]),]
cat(nrow(res2),"\n")
r<-rbind(r,res2)}

write.table(r,"probes16_aug19.txt",sep="\t",col.names=T,row.names=F,quote=F)

snps<-fread("/newshared/godmc/database_files/snps.csv",he=T)

snps2<-snps[which(snps$name%in%r$snp),]
length(snps2$name)
length(unique(snps2$name))

m<-match(r$snp,snps2$name)
df<-data.frame(r,snps2[m,"rsid"])
w<-which(is.na(df$rsid))
unique(df[w,"snp"])

write.table(df,"probes16_aug19.txt",sep="\t",quote=F,row.names=F,col.names=T)