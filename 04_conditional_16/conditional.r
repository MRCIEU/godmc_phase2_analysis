library(tidyverse)
do_conditional <- function(pval_file, bfile, pval_threshold)
{

	cmd <- paste0(
		"./gcta64 --bfile ", bfile, 
		" --cojo-file ", pval_file,
		" --cojo-joint ",
		" --cojo-p ", pval_threshold,
		" --out ", pval_file
	)
	system(cmd)
	fn <- paste0(pval_file, ".jma.cojo")
	if(!file.exists(fn)) return(NULL)

	res <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
	return(res$SNP)
}


##

arguments <- commandArgs(T)
i <- as.numeric(arguments[1])
out <- arguments[2]

bfile <- "/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur"

##

# Get cpg positions
load("cpg_pos.rdata")

a <- read_tsv(paste0("../results/16/16_", i, ".txt.gz"))
a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
a$snp2 <- a$snp
a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
a$snppos <- as.numeric(a$snppos)
a <- inner_join(a, cpgpos, by=c("cpg"))
a$cis <- FALSE
cis_radius <- 1000000
a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE

snplist <- unique(a$snp)
newbfile <- paste0("../scratch/refc_", i)
write.table(snplist, file=paste0(newbfile, ".snplist"), row=FALSE, col=FALSE, qu=FALSE)
cmd <- paste0("plink",
	" --bfile ", bfile,
	" --extract ", paste0(newbfile, ".snplist"),
	" --make-bed ", 
	" --out ", newbfile
)
system(cmd)

# Split into cis and trans
# Apply different thresholds to cis and trans
# e.g. 1e-8 for cis
# e.g. 1e-14 for trans

clumped <- group_by(a, cpg, cis) %>%
	do({
		x <- .

		# write clump file
		cpgname <- x$cpg[1]
		fn <- paste0(cpgname, "_", i, ".txt")
		system(paste0("rm ", fn, "*"))

		y <- x[, c("snp", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value")]
		y$N <- 2000
		write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)

		# Get cis/trans clumping threshold
		thresh <- ifelse(x$cis[1], 1e-4, 5e-8)
		keep <- do_conditional(fn, newbfile, thresh)
		system(paste0("rm ", fn, "*"))

		x <- subset(x, snp %in% keep)
		return(x)
	})

save(clumped, file=out)

system(paste0("rm ", newbfile, ".*"))
