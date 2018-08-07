# Library
library(data.table)
library(tidyverse)

setwd("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/nonoverlaps")
 
#create datafram of all overlaps
bait_data <-setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "CpG", "SNP", "code"))
oe_data <-setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "CpG", "SNP", "code"))

 
file_list <- list.files()
 
for (file in file_list){
    load(file)
    df1 <-data[[1]]
    df2 <-data[[2]]
    bait_data <-rbind(bait_data, df1)	
	oe_data <-rbind(oe_data, df2)	
    write.table(bait_data, file="bait_data.tsv", sep="\t", row.names=F, quote=F)
    write.table(oe_data, file="oe_data.tsv", sep="\t", row.names=F, quote=F)
}

#remove duplicates
oe_data <- subset(oe_data, select=c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "CpG", "SNP", "code"))
all_data <-rbind(bait_data, oe_data)
write.table(all_data, file="all_data.tsv", sep="\t", row.names=F, quote=F)
save(all_data, file="all.data.Rdata") # 1176 contacts overlaping 

nodups_data <- subset(all_data, !duplicated(code)) # 638 unique mQTLs (interchrom)
write.table(nodups_data, file="nodups_data.tsv", sep="\t", row.names=F, quote=F)
save(nodups_data, file="nodups.data.Rdata") 
