make_gwama <- function(out_file)
{
	message("Reading from ", out_file)
	zz <- gzfile(out_file, "rb")
	sample_size <- readBin(zz, integer(), 1)
	message("Sample size: ", sample_size)
	ncpg <- readBin(zz, integer(), 1)
	message("Number of CpGs: ", ncpg)
	cpgid <- readBin(zz, character(), ncpg)
	nsnp <- readBin(zz, integer(), 1)
	message("Number of SNPs: ", nsnp)
	snpid <- readBin(zz, character(), nsnp)
	a1 <- readBin(zz, character(), nsnp)
	a2 <- readBin(zz, character(), nsnp)
	maf <- readBin(zz, numeric(), nsnp)
	message("Expected number of associations: ", nsnp * ncpg)
	nres <- readBin(zz, integer(), 1)
	message("Reading number of associations: ", nres)
	beta <- readBin(zz, numeric(), nres, size=4)
	se <- readBin(zz, numeric(), nres, size=4)
	close(zz)

	message("Combining data")
	dat <- expand.grid(cpgid=cpgid, MARKERNAME=snpid)
	message("adding data")
	index <- match(dat$MARKERNAME, snpid)
	dat$EA <- a1[index]
	dat$NEA <- a2[index]
	dat$EAF <- maf[index]
	dat$BETA <- beta
	dat$SE <- se
	dat$N <- sample_size
	message("making markername")
	dat$MARKERNAME <- paste0(dat$MARKERNAME, "_", dat$cpgid)
	message("removing cpg col")
	dat <- subset(dat, select=-c(cpgid))
	return(dat)
}

arguments <- commandArgs(T)

fn <- arguments[1]
out <- arguments[2]

dat <- try(make_gwama(fn))
if(!class(dat) == "try-error")
{
	gz1 <- gzfile(out, "w")
	write.table(dat, file=gz1, row=FALSE, col=TRUE, qu=FALSE)
	close(gz1)
}
