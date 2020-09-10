library(data.table)
library(simulateGP)

get_ld <- function(chr, from, to, bfile, plink_bin=genetics.binaRies::get_plink_binary())
{
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()

	fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --chr ", chr,
		" --from-bp ", from, 
		" --to-bp ", to,
		" --r square ", 
		" --make-just-bim ",
		" --freq ",
		" --out ", shQuote(fn, type=shell)
	)
	system(fun1)

	x <- data.table::fread(paste0(fn, ".ld")) %>% as.matrix()
	y <- data.table::fread(paste0(fn, ".bim")) %>% as_tibble()
	z <- data.table::fread(paste0(fn, ".frq")) %>% as_tibble()
	names(y) <- c("chr", "pos", "gp", "bp", "a1", "a2")
	y$freq <- z$MAF
	unlink(paste0(fn, c(".ld", ".bim", ".frq")))
	return(list(ld=x, map=y))
}

ldetect <- "../../data/ldetect/EUR.bed"
reference <- "../../data/ref/out_hrc"

set.seed(1234)
a <- fread(ldetect) %>%
	slice(sample(1:nrow(.), 2, replace=FALSE)) %>%
	mutate(chr = gsub("chr", "", chr))
a

ld <- lapply(1:nrow(a), function(i)
{
	get_ld(a$chr[i], a$start[i], a$stop[i], reference)
})

save(ld, file="../data/ld.rdata")
