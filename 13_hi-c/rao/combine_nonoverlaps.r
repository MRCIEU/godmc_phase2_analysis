# Library
library(data.table)
library(tidyverse)

arguments <- commandArgs(T)
pnum <- (arguments[1])

print(pnum)

setwd(paste0("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/nonoverlaps/perm",pnum))
 
#create datafram of all overlaps
bait_data <-setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "cpg_samp", "SNP", "code"))
oe_data <-setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "cpg_samp", "SNP", "code"))

 
file_list <- list.files()
 
for (file in file_list){
    load(file)
    df1 <-data[[1]]
    df2 <-data[[2]]
    bait_data <-rbind(bait_data, df1)	
	oe_data <-rbind(oe_data, df2)	
    write.table(bait_data, file=paste0("bait_data_perm_",pnum,".tsv"), sep="\t", row.names=F, quote=F)
    write.table(oe_data, file=paste0("oe_data_perm_",pnum,".tsv"), sep="\t", row.names=F, quote=F)
}

#remove duplicates
oe_data <- subset(oe_data, select=c("interaction", "bait.start", "bait.end", "oe.start", "oe.end", "contacts", "chr.bait", "chr.oe", "cpg_samp", "SNP", "code"))
all_data <-rbind(bait_data, oe_data)
write.table(all_data, file=paste0("all_data_perm_",pnum,".tsv"), sep="\t", row.names=F, quote=F)
save(all_data, file=paste0("all_data_perm_",pnum,".Rdata"))

nodups_data <- subset(all_data, !duplicated(code))
write.table(nodups_data, file=paste0("nodups_data_perm_",pnum,".tsv"), sep="\t", row.names=F, quote=F)
save(nodups_data, file=paste0("nodups_data_perm_",pnum,".Rdata")) 
