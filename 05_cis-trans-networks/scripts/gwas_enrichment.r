library(dplyr)
library(data.table)
library(parallel)

perms <- function(target_pvals, background_pvals)
{
	y <- c(rep(1,length(target_pvals)), rep(0, length(background_pvals)))
	x <- -log10(c(target_pvals, background_pvals))
	x[is.infinite(x)] <- max(x[is.finite(x)])
	o <- summary(glm(y ~ x, family="binomial"))
	return(as_data_frame(coefficients(o)[2,,drop=FALSE]))
}

run <- function(j)
{
	message(j, " of ", nrow(dat))
	snps <- unique(subset(entities, cluster == dat$clust[j])$snp_rsid)
	target_pvals <- na.omit(subset(extracted, mrbaseid == dat$id[j] & ID %in% snps)$PVAL)
	background_pvals <- na.omit(subset(extracted, mrbaseid == dat$id[j] & ! ID %in% snps)$PVAL)
	if(length(target_pvals) < 2)
	{
		return(NULL)
	} else {
		o <- perms(target_pvals, background_pvals)
		o$nsnp <- length(target_pvals)
		o$j <- j
		return(o)
	}	
}

load("../data/extract_gwas.rdata")
load("../data/entity_info.rdata")

temp1 <- subset(entities, type != "creg_cpg" & snp_rsid %in% extracted$ID)
temp2 <- temp1 %>%
	group_by(cluster) %>%
	summarise(n=n()) %>%
	filter(n >= 10)

dat <- expand.grid(
	clust=unique(temp2$cluster), 
	id=unique(extracted$mrbaseid),
	nsnp=NA,
	min_p=NA,
	fisher=NA
)

o <- mclapply(1:nrow(dat), run, mc.cores=20)
save(o, entities, file="../results/gwas_enrichment.rdata")



