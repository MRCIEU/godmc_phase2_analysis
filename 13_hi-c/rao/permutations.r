# Permutation datasets

# Library
library(GenomicRanges)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationHub)
library(data.table)
library(tidyverse)

#1. Load and format data clumped data
#2. Select GRanges for SNP+proxies
#3. Gen CpG positions
#4. Sample CpGs and join correct positions
#5. remove duplicated codes
#6. loop for 4 & 5, select GRanges for CpGs and save datasets

#1.
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE )) #271724, 23117 trans, 248607 cis
clumped$code <- paste(clumped$cpg, clumped$snp)
trans_mqtl <- subset(clumped, clumped$cpgchr != clumped$snpchr) # 18584 trans mQTLs (interchromosomal)


#2.
load("/panfs/panasas01/shared-godmc/1kg_reference_ph3/snpcontrolsets_selection.rdata") #284819, 24443 trans, 260376 cis
ldinfo <- subset(f.all, select=c(SNP, min, max, nproxies)) #10085072 obs
names(ldinfo) <- c("snp", "snplow", "snphigh", "snpproxies")
temp <- subset(trans_mqtl, !duplicated(snp), select=c("snp", "snpchr", "snppos", "cpg", "cpgchr", "cpgpos", "code")) #12252 SNPs (unique)
temp <- inner_join(temp, ldinfo, "snp") #12252 obs
snps <- GRanges(seqnames=temp$snpchr, ranges=IRanges(temp$snplow, temp$snphigh), strand="*")
names(snps) <- temp$snp # 12252 unique SNPs
mcols(snps) <- temp

#3. 
cpg_pos <- subset(trans_mqtl, select=c("cpg", "cpgpos", "cpgchr"))
cpg_pos$cpg2 <- cpg_pos$cpg

#4. 

#n_perm <- 10
#set.seed(1234)

#temp <- subset(trans_mqtl, select=c("cpg", "cpgchr", "cpgpos", "code"))
#temp$cpg2 <- sample(temp$cpg)
#temp2 <- inner_join(temp, cpg_pos, "cpg2")

## Subset the samples to join with the SNP data
#temp3 <- subset(temp2, select=c("code", "cpg.y", "cpgpos.y", "cpgchr.y"))
#names(temp3) <- c("code", "cpg_samp", "cpg_pos_samp", "cpgchr_samp")

#perm1 <- inner_join(trans_mqtl, temp3, "code")
#perm1$code2 <- paste(perm1$cpg_samp, perm1$snp)


##5.
#perm1_codes <-perm1$code2 #sampled codes (25336)
#tmp <- subset(perm1_codes, !duplicated(perm1_codes)) #remove duplicated sampled codes (18583)
#codes <- trans_mqtl$code # mqtl codes (18584) (not sure why 1 more is duplicated in the perm codes than original) 
#tmp2 <- subset(codes, !duplicated(codes)) # no duplicates
#tmp3 <- subset(tmp2, !(tmp %in% tmp2)) # have a list of codes in tmp2 (mqtl codes) that are also in the sampled codes (n=7)

## remove replicated perm SNP-CpGs in real data
#perm2 <- subset(perm1, perm1$code %in% tmp3)

#temp <- subset(perm2, !duplicated(cpg_samp), select=c("cpg_samp", "cpgchr_samp", "cpg_pos_samp", "code", "code2")) # maybe don't duplicate it
#cpgs <- GRanges(seqnames=temp$cpgchr_samp, ranges=IRanges(temp$cpg_pos_samp, temp$cpg_pos_samp), strand="*")
#names(cpgs) <- temp$cpg_samp # There are 15756 unique CpGs
#mcols(cpgs) <- temp

#6. Loop

n_perm <- 1000
set.seed(1234)

for(i in 1:n_perm)
{
message(i)
temp <- subset(trans_mqtl, select=c("cpg", "cpgchr", "cpgpos", "code"))
temp$cpg2 <- sample(temp$cpg)
temp2 <- inner_join(temp, cpg_pos, "cpg2")
temp3 <- subset(temp2, select=c("code", "cpg.y", "cpgpos.y", "cpgchr.y"))
names(temp3) <- c("code", "cpg_samp", "cpg_pos_samp", "cpgchr_samp")
perm1 <- inner_join(trans_mqtl, temp3, "code")
perm1$code2 <- paste(perm1$cpg_samp, perm1$snp)
perm1_codes <-perm1$code2
tmp <- subset(perm1_codes, !duplicated(perm1_codes))
codes <- trans_mqtl$code
tmp2 <- subset(codes, !duplicated(codes))
tmp3 <- subset(tmp2, !(tmp %in% tmp2))
perm2 <- subset(perm1, perm1$code %in% tmp3)
save(perm2, file=paste0("data/permutations/permutations", i, ".rdata"))
}

