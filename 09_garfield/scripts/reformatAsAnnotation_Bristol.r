## read in meta-analysis results and reformat for use as annotation in garfield
library(tidyverse)
library(data.table)
setwd("~/repo/godmc_phase2_analysis") ## set to location where meta analysis results are

log_con <- file("~/repo/godmc_phase2_analysis/09_garfield/CisTransThresAnno.log")


uk10k<-NULL
for(chr in 1:22){
	tmp<-fread(paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/maftssd/chr", chr, sep = ""))
	tmp<-cbind(chr, tmp$V1)
	uk10k<-rbind(uk10k, tmp)
	}
cat(paste("Number of variants loaded from UK10K background files:", nrow(uk10k)), file = log_con)
	
colnames(uk10k)<-c("snpchr", "snppos")
uk10k<-as.data.frame(uk10k)
uk10k.merge<-paste("chr", uk10k[,1], ":", uk10k[,2], sep = "")
	
P_threshold<-c(1e-6,1e-8,1e-10,1e-12,1e-14)
cis_radius<-1000000
a.minP<-NULL
for(i in 1:962){
cat(i,"\n")
	load(paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_", i, ".rdata"))
	
	## filter to those in uk10k reference
	res<-mutate(res, chrPos = paste(snpchr, snppos, sep = ":"))
	res<-filter(res, chrPos %in% uk10k.merge)

	## for each batch find the minimum p for each variant
	a.minP<-rbind(a.minP, summarize(group_by(res, snpchr, snppos, cis), minP=min(Pvalue)))
	
	## as variants likely to be located across batches repeat after merge
	a.minP<-summarize(group_by(a.minP, snpchr, snppos, cis), minP=min(minP))
}

a.minP.cis<-filter(a.minP, cis == TRUE)
a.minP.trans<-filter(a.minP, cis == FALSE)

cat(paste("Number of cis variants loaded from QTL meta-analysis :", nrow(a.minP.cis)), file = log_con)
cat(paste("Number of trans variants loaded from QTL meta-analysis :", nrow(a.minP.trans)), file = log_con)
cat(paste("Number of cis & trans variants loaded from QTL meta-analysis :", nrow(intersect(select(a.minP.cis, snpchr,snppos), select(a.minP.trans, snpchr,snppos)))), file = log_con)


a.minP.cis<-mutate(a.minP.cis, annoColsCis = paste(rep(0, length(P_threshold)), collapse = ""))
a.minP.trans<-mutate(a.minP.trans, annoColsTrans = paste(rep(0, length(P_threshold)), collapse = ""))
for(j in 1:length(P_threshold)){
	a.minP.cis[which(a.minP.cis$minP < P_threshold[j]) ,"annoColsCis"]<-paste(c(rep(1, j), rep(0,(length(P_threshold)-j))),  collapse = "")
	a.minP.trans[which(a.minP.trans$minP < P_threshold[j]) ,"annoColsTrans"]<-paste(c(rep(1, j), rep(0,(length(P_threshold)-j))),  collapse = "")
}


a.minP.cis<-select(a.minP.cis, snpchr, snppos, annoColsCis)
a.minP.trans<-select(a.minP.trans, snpchr, snppos, annoColsTrans)

## add in uk10k variants not tested
for(chr in 1:22){
	uk10k.sub<-filter(uk10k, snpchr == chr)
	
	##cis
	a.sub<-filter(a.minP.cis, snpchr == paste("chr", chr, sep = ""))
	uk10k.sub<-uk10k.sub[!uk10k.sub$snppos %in% a.sub$snppos,]
	uk10k.sub<-mutate(uk10k.sub, annoColsCis = "00000")
	
	a.sub<-rbind(data.frame(select(ungroup(a.sub), snppos, annoColsCis)), data.frame(select(uk10k.sub,snppos, annoColsCis)))
	a.sub.cis<-a.sub[order(a.sub$snppos),]
	
	##trans
	a.sub<-filter(a.minP.trans, snpchr == paste("chr", chr, sep = ""))
	uk10k.sub<-uk10k.sub[!uk10k.sub$snppos %in% a.sub$snppos,]
	uk10k.sub<-mutate(uk10k.sub, annoColsTrans = "00000")
	
	a.sub<-rbind(data.frame(select(ungroup(a.sub), snppos, annoColsTrans)), data.frame(select(uk10k.sub,snppos, annoColsTrans)))
	a.sub.trans<-a.sub[order(a.sub$snppos),]

	a.minP.comb<-full_join(a.sub.cis, a.sub.trans)

	a.minP.comb<-ungroup(a.minP.comb)
	a.minP.comb$annoColsCis[which(is.na(a.minP.comb$annoColsCis))]<-paste(rep(0, length(P_threshold)), collapse = "")
	a.minP.comb$annoColsTrans[which(is.na(a.minP.comb$annoColsTrans))]<-paste(rep(0, length(P_threshold)), collapse = "")

	a.minP.comb<-mutate(a.minP.comb,
		annoCols = paste0(annoColsCis, annoColsTrans,  sep=""))
	
	
	write.table(a.minP.comb[,c(1,4)],paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/09_garfield/annotationCisTrans/chr", chr, sep = ""), col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
}


## create link file

linkFile<-cbind(0:(2*length(P_threshold)-1), c(paste("cis", P_threshold, sep = "_"), paste("trans", P_threshold, sep = "_")), "NA", "NA", c(rep("cis", length(P_threshold)), rep("trans", length(P_threshold))), "mQTL")

colnames(linkFile)<-c("Index","Annotation","Celltype","Tissue","Type","Category")

write.table(linkFile, "/panfs/panasas01/shared-godmc/godmc_phase2_analysis/09_garfield/annotationCisTrans/link_file.txt", row.names = FALSE, sep = " ", quote = FALSE)

