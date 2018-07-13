# 1. Filter communities so that there are no correlated CpGs due to physical proximity

library(LOLA)
library(dplyr)
library(GenomicRanges)
library(doParallel)
library(rvest)

load("../data/trans_granges.rdata")

# Read in stuff
tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

anno <- tfbsdb$regionAnno

# Register parallel
# (no_cores <- detectCores() - 1)

no_cores <- 10
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores, type="FORK")

# Extract SNPs
snpres <- parLapply(cl, 1:nrow(anno), function(i)
{
	message(i, " of ", nrow(anno))
	temp <- names(snp)[queryHits(findOverlaps(snp, tfbsdb$regionGRL[[i]]))]
	if(length(temp) > 0)
		return(tibble(anno=i, snp=temp))
	else
		return(tibble())
})

# Extract CpGs
cpgres <- parLapply(cl, 1:nrow(anno), function(i)
{
	message(i, " of ", nrow(anno))
	temp <- names(cpg)[queryHits(findOverlaps(cpg, tfbsdb$regionGRL[[i]]))]
	if(length(temp) > 0)
		return(tibble(anno=i, cpg=temp))
	else
		return(tibble())
})


stopCluster(cl)
snpres <- bind_rows(snpres)
cpgres <- bind_rows(cpgres)

save(snpres, cpgres, anno, file="../data/annotations.rdata")



# Get tissue types for encode

t1 <- read_html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=1&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t2 <- read_html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=2&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t3 <- read_html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=3&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()

encode <- bind_rows(t1, t2, t3)

# Index annotations
anno$index <- 1:nrow(anno)

# Everything without a tissue but with a celltype
temp1 <- subset(anno, is.na(tissue) & !is.na(cellType))
encodet <- subset(encode, select=c(cell, Tissue))
temp1$tissue <- encodet$Tissue[
	match(gsub("-", "", tolower(temp1$cell)), gsub("-", "", tolower(encodet$cell)))
]

table(is.na(temp1$tissue), temp1$collection)

# Merge back

temp2 <- subset(anno, !index %in% temp1$index)
anno2 <- rbind(temp1, temp2)
anno2 <- anno2[order(anno2$index),]
stopifnot(all(anno2$filename == anno$filename))
anno <- anno2

anno$antibody[is.na(anno$antibody) & anno$collection == "sheffield_dnase"] <- anno$description[is.na(anno$antibody) & anno$collection == "sheffield_dnase"]

anno$antibody[is.na(anno$antibody) & anno$collection == "encode_segmentation"] <- anno$description[is.na(anno$antibody) & anno$collection == "encode_segmentation"]

anno$antibody[is.na(anno$antibody) & anno$collection == "ucsc_features"] <- anno$description[is.na(anno$antibody) & anno$collection == "ucsc_features"]


anno$antibody[anno$collection == "sheffield_dnase"] <- NA

anno$antibody2 <- sapply(strsplit(as.character(anno$antibody), split="_"), function(x) x[1])
anno$antibody2 <- sapply(strsplit(as.character(anno$antibody2), split=" "), function(x) x[1])
anno$antibody2 <- sapply(strsplit(as.character(anno$antibody2), split="\\("), function(x) x[1])
anno$antibody2 <- gsub("eGFP-", "", anno$antibody2)


difres$snpanno <- anno$antibody2[difres$Var1]
difres$cpganno <- anno$antibody2[difres$Var2]
difres$sddif2 <- difres$sddif
difres$sddif2[difres$val < difres$Mean] <- difres$sddif2[difres$val < difres$Mean] * -1

summary(difres$sddif2)
# hist(difres$sddif2)
# All values that are 'depleted' are not significant

blood <- subset(anno, tissue == "blood")

save(blood, anno, file="../data/blood.rdata")
