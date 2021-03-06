library(dplyr)
library(tidyr)
library(data.table)
do_conditional <- function(pval_file, bfile, pval_threshold)
{

	cmd <- paste0(
		"./gcta64 --bfile ", bfile, 
		" --cojo-file ", pval_file,
		" --cojo-slct ",
		" --cojo-p ", pval_threshold,
		" --out ", pval_file
	)
	system(cmd)
	fn <- paste0(pval_file, ".jma.cojo")
	if(!file.exists(fn)) return(data.frame())

	res <- read.table(fn, header=TRUE,
		colClasses=c("integer", "character", "integer", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")) %>% as.tbl

	# res <- read_tsv(fn,
	# 	col_types=cols(
	# 	  Chr = col_integer(),
	# 	  SNP = col_character(),
	# 	  bp = col_integer(),
	# 	  refA = col_character(),
	# 	  freq = col_double(),
	# 	  b = col_double(),
	# 	  se = col_double(),
	# 	  p = col_double(),
	# 	  n = col_double(),
	# 	  freq_geno = col_double(),
	# 	  bJ = col_double(),
	# 	  bJ_se = col_double(),
	# 	  pJ = col_double(),
	# 	  LD_r = col_double()
	# 	)
	# )
	return(res)
}


# group_by(a, Var1) %>% 
# do({ 
# 	x <- .
# 	if(x$Var1 == 3)
# 		return(data.frame())
# 	else
# 		return(as.tbl(x))
# })

##
setwd("/mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/04_conditional_16")

arguments <- commandArgs(T)
i <- as.numeric(arguments[1])
out.cis <- arguments[2]
out.trans <- arguments[3]

bfile <- "../data/ref/out_hrc"

##

# Get cpg positions
load(paste0("../results_3coh/16_3cohorts/16_cleaned_", i, ".rdata"))


snplist <- unique(res$snp)
newbfile <- paste0("../scratch_3coh/refc_", i)
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
# res <- head(res, 100)

#clumped <- group_by(res, cpg, cis) %>%
#	do({
#		x <- .
		# write clump file
#		cpgname <- x$cpg[1]
#		fn <- paste0(cpgname, "_", i, ".txt")
#		system(paste0("rm -f ", fn, "*"))

#		y <- x[, c("snp", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize")]
#		write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)

cis<-res[which(res$cis=="TRUE"),]
clumped<-data.frame()
cpgs<-unique(res$cpg)

for (j in 1:length(cpgs)){
cpgname <- cpgs[j]
fn <- paste0(cpgname, "_", i, ".txt")
system(paste0("rm -f ", fn, "*"))	
y<-cis[cis$cpg==cpgs[j],c("snp", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize")]
write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)



		# Get cis/trans clumping threshold
		#thresh <- ifelse(x$cis[1], 1e-4, 5e-8)
		thresh <-1e-4
		cat(fn,"\n")
		keep <- do_conditional(fn, newbfile, thresh)
		system(paste0(
			"cp ", fn, "*bad* bad_3coh/"
		))
        
		system(paste0("rm ", fn, "*"))
        
        keep<-data.frame(cpg=rep(cpgs[j],nrow(keep)),keep)
		clumped<-rbind(clumped,keep)
		#return(keep)
	}
#system(paste0("rm ", newbfile, ".*"))
save(clumped, file=out.cis)

####
trans<-res[which(res$cis=="FALSE"),]
clumped<-data.frame()
cpgs<-unique(res$cpg)

for (j in 1:length(cpgs)){
cpgname <- cpgs[j]
fn <- paste0(cpgname, "_", i, ".txt")
system(paste0("rm -f ", fn, "*"))	
y<-cis[cis$cpg==cpgs[j],c("snp", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize")]
write.table(y, file=fn, row=FALSE, col=TRUE, qu=FALSE)



		# Get cis/trans clumping threshold
		#thresh <- ifelse(x$cis[1], 1e-4, 5e-8)
		thresh <-5e-8
		cat(fn,"\n")
		keep <- do_conditional(fn, newbfile, thresh)
		system(paste0(
			"cp ", fn, "*bad* bad_3coh/"
		))
        
		system(paste0("rm ", fn, "*"))
        
        keep<-data.frame(cpg=rep(cpgs[j],nrow(keep)),keep)
		clumped<-rbind(clumped,keep)
		#return(keep)
	}
system(paste0("rm ", newbfile, ".*"))
save(clumped, file=out.trans)






