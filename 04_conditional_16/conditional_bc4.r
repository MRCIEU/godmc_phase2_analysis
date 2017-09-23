library(tidyverse)
do_conditional <- function(pval_file, bfile, pval_threshold)
{

	cmd <- paste0(
		"./gcta64 --bfile ", bfile, 
		" --cojo-file ", pval_file,
		" --cojo-slct ",
		" --cojo-p ", pval_threshold,
		" --out ", pval_file
	)
	system(cmd, show.output.on.console = FALSE)
	fn <- paste0(pval_file, ".jma.cojo")
	if(!file.exists(fn)) return(NULL)

	res <- read_tsv(fn,
		col_types=cols(
		  Chr = col_integer(),
		  SNP = col_character(),
		  bp = col_integer(),
		  refA = col_character(),
		  freq = col_double(),
		  b = col_double(),
		  se = col_double(),
		  p = col_double(),
		  n = col_double(),
		  freq_geno = col_double(),
		  bJ = col_double(),
		  bJ_se = col_double(),
		  pJ = col_double(),
		  LD_r = col_double()
		)
	)
	return(res)
}


##

arguments <- commandArgs(T)
i <- as.numeric(arguments[1])
out <- arguments[2]

bfile <- "../data/ref/eur"

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
a <- subset(a, (cis & Pvalue < 1e-4) | (!cis & Pvalue < 5e-8))

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
		system(paste0("rm -f ", fn, "*"))

		y <- x[, c("snp", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize")]
		write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)

		# Get cis/trans clumping threshold
		thresh <- ifelse(x$cis[1], 1e-4, 5e-8)
		keep <- do_conditional(fn, newbfile, thresh)
		print(head(keep))
		system(paste0("rm ", fn, "*"))
		return(keep)
	})

save(clumped, file=out)

system(paste0("rm ", newbfile, ".*"))
