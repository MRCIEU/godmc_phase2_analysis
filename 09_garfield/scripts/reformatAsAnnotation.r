## read in meta-analysis results and reformat for use as annotation in garfield
library(tidyverse)
setwd("/mnt/data1/goDMC_Phase2/") ## set to location wher meta16 is


P_threshold<-c(5e-5,5e-6,5e-7,5e-8,5e-9,5e-10,5e-11)
cis_radius<-500000
a.minP<-NULL
for(i in 1:962){

	a <- read_tsv(paste0("meta16/16_", i, ".txt.gz"))
	a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
	a$snp2 <- a$snp
	a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
	a$snppos <- as.numeric(a$snppos)
	#a <- inner_join(a, cpgpos, by=c("cpg"))
	#a$cis <- FALSE
	#a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE
	## for each batch find the minimum p for each variant
	a.minP<-rbind(a.minP, summarize(group_by(a, snpchr, snppos), minP=min(Pvalue)))
	
	## as variants likely to be located across batches repeat after merge
	a.minP<-summarize(group_by(a.minP, snpchr, snppos), minP=min(minP))
}

a.minP<-mutate(a.minP, annoCols = paste(rep(0, length(P_threshold)), collapse = ""))
for(j in 1:length(P_threshold)){
	a.minP[which(a.minP$minP < P_threshold[j]) ,"annoCols"]<-paste(c(rep(1, j), rep(0,(length(P_threshold)-j))),  collapse = "")
}

write.csv(table(a.minP$annoCols), "godmc_phase2_analysis/09_garfield/results/TableSNPsAnnotations.csv")

a.minP<-ungroup(a.minP)
for(chr in 1:22){
	a.sub<-filter(a.minP, snpchr == paste("chr", chr, sep = ""))
	
	write.table(select(a.sub,snppos, annoCols),paste("godmc_phase2_analysis/09_garfield/annotation/chr", chr, sep = ""), col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
}
