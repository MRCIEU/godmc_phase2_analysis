suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

flip_alleles <- function(dat, fliplist)
{
	index <- dat$SNP %in% fliplist
	EA <- dat$EA
	NEA <- dat$NEA
	EA[index & dat$EA == "A"] <- "T"
	EA[index & dat$EA == "C"] <- "G"
	EA[index & dat$EA == "T"] <- "A"
	EA[index & dat$EA == "G"] <- "C"
	NEA[index & dat$NEA == "A"] <- "T"
	NEA[index & dat$NEA == "C"] <- "G"
	NEA[index & dat$NEA == "T"] <- "A"
	NEA[index & dat$NEA == "G"] <- "C"
	dat$EA <- EA
	dat$NEA <- NEA
	return(dat)
}

check_alleles <- function(dat, bim)
{
	temp <- subset(bim, snp %in% dat$SNP)
	temp <- merge(dat, temp, by.x="SNP", by.y="snp")
	return(table(
		(temp$EA == temp$a1 & temp$NEA == temp$a2) |
		(temp$EA == temp$a2 & temp$NEA == temp$a1)
	))
}


fn <- commandArgs(T)[1]
flips <- commandArgs(T)[2]
ref <- commandArgs(T)[3]

# load(ref)
# bim$a1 <- toupper(bim$a1)
# bim$a2 <- toupper(bim$a2)

a <- fread(paste0("zcat ", fn))
fl <- scan(flips, what="character")

# check_alleles(a, bim)
dat <- flip_alleles(a, fl)
# check_alleles(dat, bim)

message("writing file to ", fn)
con <- gzfile(fn, "w")
write.table(dat, con, row=F, col=T, qu=F)
close(con)
message("done")

# b <- subset(a, SNP %in% fl)
# bim2 <- subset(bim, snp %in% b$SNP)

# temp <- merge(b, bim2, by.x="SNP", by.y="snp")

# temp$a1 <- toupper(temp$a1)
# temp$a2 <- toupper(temp$a2)

# table(temp$)
