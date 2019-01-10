library(TwoSampleMR)
library(MRInstruments)
library(dplyr)

load("../../05_cis-trans-networks/data/labelids.rdata")
data(mrbase_instruments)

gwas <- subset(mrbase_instruments, id.exposure %in% labelids$id)
ids <- unique(gwas$id.exposure)
dat <- data_frame(id=ids, maxr2=NA, meanr2=NA)

for(i in 1:length(ids))
{
	message(i)
	temp <- subset(gwas, id.exposure == ids[i])$SNP
	write.table(temp, file="temp", row=F, col=F, qu=F)
	cmd <- paste0(
		"plink --bfile ~/repo/mr-base-api/app/ld_files/data_maf0.01_rs --r2 square --extract temp --out temp"
	)
	system(cmd)
	a <- read.table("temp.ld",header=FALSE) %>% as.matrix
	a <- a[lower.tri(a)]
	dat$maxr2[i] <- max(a)
	dat$meanr2[i] <- mean(a)
}

hist(dat$maxr2)
hist(dat$meanr2)

save(dat, file="../data/instrument_r2.rdata")

