#!/bin/bash


for i in {1..300}
do
echo $i
zfgrep -f ~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/metab.txt ../../results/17/old/17_${i}.txt.gz > temp${i}
done

i=10

parallel ./westra.sh ::: {1..300}
parallel ./metab.sh ::: {73..300}

fgrep -wf ~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/gtex.txt ../../../../1kg_reference_ph3/eur.filtered.bim | wc -l

parallel ./gtex.sh ::: {1..300}

sewlibrary(data.table)
library(tidyverse)


l <- list()
for(i in 1:300)
{
	a <- try(fread(paste0("westra", i)))
	if(class(a) != 'try-error')
	{	
		if(nrow(a) > 0)
		{
			message(i)
			a <- tidyr::separate(a, V1, c("snp", "cpg"), sep="_")
			a <- subset(a, V8 < 1e-3)
			l[[i]] <- a
		}
	}
}



library(data.table)
lg <- list()
for(i in 1:300)
{
	a <- try(fread(paste0("temp", i)))
	if(class(a) != 'try-error')
	{	
		if(nrow(a) > 0)
		{
			message(i)
			a <- tidyr::separate(a, V1, c("snp", "cpg"), sep="_")
			a <- subset(a, V8 < 1e-3)
			lg[[i]] <- a
		}
	}
}

library(tidyverse)

for(i in 1:length(lg))
{
	if(!is.null(lg[[i]]))
		{
	if(nrow(lg[[i]]) == 0)
		lg[[i]] <- NULL

		}
}

lg2 <- plyr::rbind.fill(lg)

lg3 <- subset(lg2, V8 < 1.2e-10)

dim(lg)

library(meffil)

pi <- meffil.get.features() %>% dplyr::select(name, chromosome, position)
lg <- inner_join(lg, pi, by=c("cpg"="name"))

lg$SNP <- lg$snp
lg <- tidyr::separate(lg, snp, c("snpchr", "snppos", "snptype"), ":")
table(lg$chromosome == lg$snpchr)

tapply(-log10(lg$V8), lg$chromosome == lg$snpchr, mean)

transthresh <- 0.05 / (2500 * 350000)
cisthresh <- 0.05 / (2500 * 100)

lg$snppos <- as.numeric(lg$snppos)

lg$cis <- lg$snpchr == lg$chromosome & abs(lg$snppos - lg$position) < 1000000

tapply(-log10(lg$V8), lg$cis, mean)

lg$sig <- lg$V8 < transthresh
lg$sig[lg$cis] <- lg$V8[lg$cis] < cisthresh



temp <- do.call(c, strsplit(abc$celltype, split=","))
as.data.frame(table(temp))




