library(tidyverse)

i <- 1


do_clump <- function(pval_file, bfile, pval_threshold, rsq_threshold, kb_threshold)
{
	cmd <- paste0(
		"module add apps/plink-1.90b3v;",
		"plink --bfile ", bfile, 
		"--clump ", pval_file,
		"--clump-p1 ", pval_threshold,
		"--clump-r2 ", rsq_threshold,
		"--clump-kb ", kb_threshold,
		"--out ", pval_file
	)
	res <- read.table(pval_file, header=TRUE, stringsAsFactors=FALSE)
	return(res$SNP)
}


a <- read_tsv(paste0("../results/16/16_", i, ".txt.gz"))
a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")


# Split into cis and trans

# Apply different thresholds to cis and trans

# e.g. 1e-8 for cis
# e.g. 1e-14 for trans


a <- filtergroup_by(a, cpg) %>%
	do({
		x <- .
		# write clump file
		cpgname <- x$cpg[1]
		y <- x[, c("snp", "P-value")]
		names(y) <- c("SNP", "P")
		write.table(y, row=FALSE, col=TRUE, qu=FALSE)



	})


