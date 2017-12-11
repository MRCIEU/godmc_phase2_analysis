# Script to format Grubert Hi-C data for enrichment analysis 

# Library
library(data.table)
library(tidyverse)

arguments <- commandArgs(T)
chunk <- arguments[1]
#chunk <-"001"
message(chunk)

# Load datasets
mid <- fread(paste0("zcat chunks/grubert_",chunk,".gz", sep=""), header = T) # Interactions (mid points only)
frag <- fread("zcat fragmentStartEnd.txt.gz", header = F) # Interactions (start/end)

# subset to cor of 0.4 (matches)
mid <-subset(mid, mid$cor >=0.4)

# names
colnames(frag)[which(names(frag) == "V1")] <- "chr"
colnames(frag)[which(names(frag) == "V2")] <- "start"
colnames(frag)[which(names(frag) == "V3")] <- "end"
colnames(frag)[which(names(frag) == "V4")] <- "mid"

# Merge start and end with mid points
m1 <- merge(mid, frag, by.x = c("chr", "pos_i"), by.y = c("chr", "mid"))
colnames(m1)[which(names(m1) == "start")] <- "i.start"
colnames(m1)[which(names(m1) == "end")] <- "i.end"
m1 <- merge(m1, frag, by.x = c("chr", "pos_j"), by.y = c("chr", "mid"))
colnames(m1)[which(names(m1) == "start")] <- "j.start"
colnames(m1)[which(names(m1) == "end")] <- "j.end"
HiC <- m1[order(m1$pair.id),]
HiC <- as.data.frame(HiC)

HiC$id = paste(HiC$pos_i, HiC$pos_j, sep="_")
HiC$dist <-(HiC$pos_j - HiC$pos_i)
HiC <-HiC[c("id", "chr", "i.start", "i.end", "j.start", "j.end", "pos_i", "pos_j", "dist", "i", "j", "cor", "x")]

write.table(HiC, file=gzfile(paste0("chunks/chunks_clean/grubert_cor0.4_",chunk,".gz", sep="")), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(HiC, file=paste0("chunks/chunks_clean/grubert_cor0.4_",chunk,".tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


