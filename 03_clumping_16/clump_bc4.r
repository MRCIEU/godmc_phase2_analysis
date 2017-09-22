library(dplyr)
library(data.table)
library(tidyr)
do_clump <- function(pval_file, bfile, pval_threshold, rsq_threshold, kb_threshold)
{
	cmd <- paste0(
		"plink --bfile ", bfile, 
		" --clump ", pval_file,
		" --clump-p1 ", pval_threshold,
		" --clump-r2 ", rsq_threshold,
		" --clump-kb ", kb_threshold,
		" --out ", pval_file
	)
	system(cmd)
	fn <- paste0(pval_file, ".clumped")
	if(!file.exists(fn)) return(NULL)

	res <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
	return(res$SNP)
}


##

arguments <- commandArgs(T)
i <- as.numeric(arguments[1])
out <- arguments[2]
cis_pval <- as.numeric(arguments[3])
trans_pval <- as.numeric(arguments[4])
rsq_thresh <- as.numeric(arguments[5])
kb_thresh <- as.numeric(arguments[6])
cis_radius <- as.numeric(arguments[7])

#bfile <- "/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur"
bfile <- "data/eur"
##

# Get cpg positions
load("cpg_pos.rdata")

a <- fread(paste0("zcat ../results/16/16_", i, ".txt.gz"))
a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
a$snp2 <- a$snp
a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
a$snppos <- as.numeric(a$snppos)
a <- inner_join(a, cpgpos, by=c("cpg"))
a$cis <- FALSE
a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE


# Make specific file for each run
# Might improve speed / reduce IO problem?
snplist <- unique(a$snp)
newbfile <- paste0("../scratch/ref_", i)
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


clumped <- group_by(a, cpg, cis) %>%
	do({
		x <- .

		# write clump file
		cpgname <- x$cpg[1]
		fn <- paste0(cpgname, "_", i, ".txt")
		system(paste0("rm -f ", fn, "*"))

		y <- x[, c("snp", "Pvalue")]
		names(y) <- c("SNP", "P")
		write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)

		# Get cis/trans clumping threshold
		thresh <- ifelse(x$cis[1], cis_pval, trans_pval)
		keep <- do_clump(fn, newbfile, thresh, rsq_thresh, kb_thresh)
		system(paste0("rm ", fn, "*"))

		x <- subset(x, snp %in% keep)
		return(x)
	})

save(clumped, file=out)

system(paste0("rm ", newbfile, ".*"))


