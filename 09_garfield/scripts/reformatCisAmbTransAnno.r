### use files categorizing snps as cis, trans or ambivalent QTL to create annotation files for garfield to do GWAS enrichment


library(tidyverse)
library(data.table)
setwd("/mnt/data1/goDMC_Phase2/") ## set to location where meta-analysis results are

## load uk10k variants for background list
uk10k<-NULL
for(chr in 1:22){
	tmp<-fread(paste("/mnt/data1/programs/garfield-data/maftssd/chr", chr, sep = ""))
	tmp<-cbind(chr, tmp$V1)
	uk10k<-rbind(uk10k, tmp)
	}
	
colnames(uk10k)<-c("snpchr", "snppos")
uk10k<-as.data.frame(uk10k)
uk10k.merge<-paste("chr", uk10k[,1], ":", uk10k[,2], sep = "")

## open snp classification files
snpFile <- read.table(gzfile("godmc_phase2_analysis/16_jan2018/snpcpgpval.chr1.cistrans.txt.gz"), header = TRUE)
for(chr in 2:22){
	snpFile<-rbind(snpFile, read.table(gzfile(paste("godmc_phase2_analysis/16_jan2018/snpcpgpval.chr", chr, ".cistrans.txt.gz", sep = "")), header = TRUE))
}

	