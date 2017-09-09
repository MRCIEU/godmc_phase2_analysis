## read in meta-analysis results and reformat for use as annotation in garfield
library(tidyverse)
setwd("/mnt/data1/goDMC_Phase2/") ## set to location wher meta16 is

load("/mnt/data1/450K_reference/AllProbeIlluminaAnno.Rdata")
probeAnnot<-probeAnnot[,c("NAME", "CHR", "MAPINFO")]
colnames(probeAnnot)<-c("cpg", "cpgchr", "cpgpos")
probeAnnot$cpg<-as.character(probeAnnot$cpg)
probeAnnot$cpgchr<-as.character(probeAnnot$cpgchr)
probeAnnot$cpgchr<-paste("chr", probeAnnot$cpgchr, sep = "")

P_threshold<-c(1e-6,1e-8,1e-10,1e-12,1e-14)
cis_radius<-1000000
a.minP<-NULL
for(i in 1:962){

	a <- read_tsv(paste0("meta16/16_", i, ".txt.gz"))
	a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
	a$snp2 <- a$snp
	a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
	a$snppos <- as.numeric(a$snppos)
	a <- inner_join(a, probeAnnot, by=c("cpg"))
	a$cis <- FALSE
	a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE
	## for each batch find the minimum p for each variant
	a.minP<-rbind(a.minP, summarize(group_by(a, snpchr, snppos, cis), minP=min(Pvalue)))
	
	## as variants likely to be located across batches repeat after merge
	a.minP<-summarize(group_by(a.minP, snpchr, snppos, cis), minP=min(minP))
}

a.minP.cis<-filter(a.minP, cis == TRUE)
a.minP.trans<-filter(a.minP, cis == FALSE)

a.minP.cis<-mutate(a.minP.cis, annoColsCis = paste(rep(0, length(P_threshold)), collapse = ""))
a.minP.trans<-mutate(a.minP.trans, annoColsTrans = paste(rep(0, length(P_threshold)), collapse = ""))
for(j in 1:length(P_threshold)){
	a.minP.cis[which(a.minP.cis$minP < P_threshold[j]) ,"annoColsCis"]<-paste(c(rep(1, j), rep(0,(length(P_threshold)-j))),  collapse = "")
	a.minP.trans[which(a.minP.trans$minP < P_threshold[j]) ,"annoColsTrans"]<-paste(c(rep(1, j), rep(0,(length(P_threshold)-j))),  collapse = "")
}

write.csv(rbind(table(a.minP.cis$annoColsCis),table(a.minP.trans$annoColsTrans)),"godmc_phase2_analysis/09_garfield/results/TableSNPsAnnotations.csv")

a.minP.cis<-select(a.minP.cis, snpchr, snppos, annoColsCis)
a.minP.trans<-select(a.minP.trans, snpchr, snppos, annoColsTrans)

## add in background SNPs that weren't tested in the meta-anlysis
bim<-fread("/mnt/data1/goDMC_Phase2/1kg_reference/eur.filtered.bim")

for(chr in paste("chr", 1:22, sep = "")){
	bp<-as.data.frame(bim[which(bim[,1] == gsub("chr", "", chr)),4])
	bp<-bp[,1]

	lookUp<-select(ungroup(filter(a.minP.cis, snpchr == chr)),  snppos)
	toAdd<-bp[!bp %in% lookUp$snppos]
	toAdd<-data.frame("snpchr" = chr, "snppos" = toAdd, "annoColsCis" = paste(rep(0, length(P_threshold)), collapse = ""))
	a.minP.cis<-full_join(a.minP.cis, toAdd)
	
	lookUp<-select(ungroup(filter(a.minP.trans, snpchr == chr)),  snppos)
	toAdd<-bp[!bp %in% lookUp$snppos]
	toAdd<-data.frame("snpchr" = chr, "snppos" = toAdd, "annoColsTrans" = paste(rep(0, length(P_threshold)), collapse = ""))
	a.minP.trans<-full_join(a.minP.trans, toAdd)
	
}

a.minP.comb<-full_join(a.minP.cis, a.minP.trans)

a.minP.comb<-ungroup(a.minP.comb)
a.minP.comb$annoColsCis[which(is.na(a.minP.comb$annoColsCis))]<-paste(rep(0, length(P_threshold)), collapse = "")
a.minP.comb$annoColsTrans[which(is.na(a.minP.comb$annoColsTrans))]<-paste(rep(0, length(P_threshold)), collapse = "")

a.minP.comb<-mutate(a.minP.comb,
   annoCols = paste0(annoColsCis, annoColsTrans,  sep=""))


for(chr in 1:22){
	a.sub<-filter(a.minP.comb, snpchr == paste("chr", chr, sep = ""))
	
	write.table(select(a.sub,snppos, annoCols),paste("godmc_phase2_analysis/09_garfield/annotation/chr", chr, sep = ""), col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
}
chr<-23
a.sub<-filter(a.minP.comb, snpchr == paste("chr", chr, sep = ""))
write.table(select(a.sub,snppos, annoCols),"godmc_phase2_analysis/09_garfield/annotation/chrX", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)


## create link file

linkFile<-cbind(0:(2*length(P_threshold)-1), c(paste("cis", P_threshold, sep = "_"), paste("trans", P_threshold, sep = "_")), "NA", "NA", c(rep("cis", length(P_threshold)), rep("trans", length(P_threshold))), "mQTL")

colnames(linkFile)<-c("Index","Annotation","Celltype","Tissue","Type","Category")

write.table(linkFile, "godmc_phase2_analysis/09_garfield/annotation/link_file.txt", row.names = FALSE, sep = " ", quote = FALSE)

