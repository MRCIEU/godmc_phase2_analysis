library(ggplot2)
library(coloc)
library(simulateGP)


ldetect <- "../data/ldetect/EUR.bed"
reference <- "../data/ref/out_hrc"


params <- expand.grid(
	n = 28000,
	region = 1,
	ncausal = c(0, 1, 2, 3),
	coverage = c(1, 0.8, 0.6, 0.4, 0.2, 0.1),
	rsq_trait1 = c(0.3, 0.1, 0.01),
	scale_trait2 = c(1, 0.5, 0.1),
	method = c("fill", "sparse"),
	nsim = c(1:100)
)












get_ld <- function(variants, bfile, plink_bin=genetics.binaRies::get_plink_binary())
{
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)

		fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --extract ", shQuote(fn, type=shell), 
		" --recode A ", 
		" --out ", shQuote(fn, type=shell)
	)
	system(fun1)

	x <- data.table::fread(paste0(fn, ".raw")) %>% {.[,-c(1:6)]} %>% as.matrix()
	unlink(paste0(fn, ".raw"))
	return(x)
}








