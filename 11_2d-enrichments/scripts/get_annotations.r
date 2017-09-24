# 1. Filter communities so that there are no correlated CpGs due to physical proximity

library(LOLA)
library(dplyr)
library(GenomicRanges)
library(doParallel)

load("../data/trans_granges.rdata")

# Read in stuff
tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")

# Register parallel
(no_cores <- detectCores() - 1)
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

save(snpres, cpgres, file="../results/annotations.rdata")
