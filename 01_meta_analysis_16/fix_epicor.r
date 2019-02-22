suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

update_alleles <- function(dat, indels)
{
	w<-which(names(dat)%in%c("A1"))
	if(length(w)==1){names(dat)[w]<-c("EA")}
	w<-which(names(dat)%in%c("A2"))
	if(length(w)==1){names(dat)[w]<-c("NEA")}

	dat$index <- 1:nrow(dat)
	da <- subset(dat, SNP %in% indels$SNP)
	da <- merge(da, indels, by="SNP")
	message("Fixing ", nrow(da), " indels")
	i1 <- da$EA== da$A1.orig
	i2 <- da$NEA == da$A2.orig
	da$EA[i1] <- da$A1[i1]
	da$NEA[i1] <- da$A2[i1]
	da$EA[i2] <- da$A1[i2]
	da$NEA[i2] <- da$A2[i2]
	da <- subset(da, select=-c(pos, A1, A2, A1.orig, A2.orig))
	dat <- subset(dat, !SNP %in% indels$SNP)
	dat <- rbind(dat, da)
	dat <- dat[order(dat$index), ]
	dat <- subset(dat, select=-c(index))
}

fn <- commandArgs(T)[1]
indels <- commandArgs(T)[2]

a <- suppressWarnings(fread(paste0("zcat ", fn)))
indels <- fread(indels)

# check_alleles(a, bim)
dat <- update_alleles(a, indels)
index <- dat$EA %in% c("Z","Y") | dat$NEA %in% c("Z", "Y")
message("Removing ", sum(index), " stragglers")
dat <- subset(dat, !index)
# check_alleles(dat, bim)

message("writing file to ", fn)
con <- gzfile(fn, "w")
write.table(dat, con, row=F, col=T, qu=F)
close(con)
message("done")
