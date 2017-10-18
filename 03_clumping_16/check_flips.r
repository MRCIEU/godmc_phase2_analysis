# Checking flipped alleles 
# The flipping happens in 02a
# Seems to be only for SNPs that have maf close to 0.5
# For non-palindromic SNPs we can find these easily by comparing to reference dataset
# For palindromic SNPs this is not possible because 
# If ref = A T with AF 0.49
#    tar = T A with AF 0.49
#    tar = A T with AF 0.51
#    tar = T A with AF 0.51
#    tar = A T with AF 0.49
# They would all look the same after flipping
# Need the list of flipped alleles to be uploaded, which can be used to infer the SNPs that have been flipped

setwd("~/sandpit/godmc/relateds")

target <- fread("input_data/data.bim")
target_adjusted <- fread("processed_data/genetic_data/data.bim")
reference <- fread("zcat processed_data/genetic_data/CLEANED.data.easyqc.gz")

freq <- fread("zcat results/02/data.frq.gz")

target$snp <- paste0("chr", target$V1, ":", target$V4, ":SNP")

target <- subset(target, snp %in% target_adjusted$V2)
table(duplicated(target$snp))
a <- target$snp[duplicated(target$snp)]
target <- subset(target, !snp %in% a)
dim(target_adjusted)
dim(target)
dim(reference)

target <- subset(target, snp %in% target_adjusted$V2)
target_adjusted <- subset(target_adjusted, V2 %in% target$snp)
reference <- subset(reference, SNP %in% target$snp)
freq <- subset(freq, SNP %in% target$snp)

ind <- match(target$snp, target_adjusted$V2)
target <- target[ind, ]
stopifnot(all(target$snp == target_adjusted$V2))

ind <- match(target$snp, reference$SNP)
reference <- reference[ind, ]
stopifnot(all(target$snp == reference$SNP))

stopifnot(all(freq$SNP == reference$SNP))

same <- target$V5 == target_adjusted$V5 & target$V6 == target_adjusted$V6
wrong_ea <- target$V5 == target_adjusted$V6 & target$V6 == target_adjusted$V5
flipped <- ! same & ! wrong_ea
table(flipped)

head(target_adjusted[flipped,])
head(target[flipped,])

table(paste(target_adjusted$V5[flipped], target_adjusted$V6[flipped]))
table(paste(target_adjusted$V5[wrong_ea], target_adjusted$V6[wrong_ea]))
table(paste(target_adjusted$V5[same], target_adjusted$V6[same]))



same <- reference$EFFECT_ALLELE == target_adjusted$V5 & reference$OTHER_ALLELE == target_adjusted$V6
wrong_ea <- reference$EFFECT_ALLELE == target_adjusted$V6 & reference$OTHER_ALLELE == target_adjusted$V5
flipped <- ! same & ! wrong_ea
table(flipped)
table(wrong_ea)

head(target_adjusted[flipped,])
head(reference[flipped,])

table(paste(target_adjusted$V5[flipped], target_adjusted$V6[flipped]))
table(paste(target_adjusted$V5[wrong_ea], target_adjusted$V6[wrong_ea]))
table(paste(target_adjusted$V5[same], target_adjusted$V6[same]))

# There are no wrong_ea because all of these have been flipped
# The non-palindromic ones are obvious to flip back
# Some of the palindromic SNPs need to be flipped back - which ones?
# Can we work this out from MAF?


same <- reference$EFFECT_ALLELE == freq$A1 & reference$OTHER_ALLELE == freq$A2
wrong_ea <- reference$EFFECT_ALLELE == freq$A2 & reference$OTHER_ALLELE == freq$A1
flipped <- ! same & ! wrong_ea
table(flipped)
table(wrong_ea)

head(target_adjusted[flipped,])
head(reference[flipped,])

table(paste(target_adjusted$V5[flipped], target_adjusted$V6[flipped]))
table(paste(target_adjusted$V5[wrong_ea], target_adjusted$V6[wrong_ea]))
table(paste(target_adjusted$V5[same], target_adjusted$V6[same]))




library(dplyr)

index <- with(clumped, 
	((Allele1 == "a" & Allele2 == "t") |
	(Allele1 == "t" & Allele2 == "a") |
	(Allele1 == "g" & Allele2 == "c") |
	(Allele1 == "c" & Allele2 == "g"))
)
table(index)

a <- subset(clumped, index & (Freq1 > 0.45 & Freq1 < 0.55))
b <- subset(clumped, !index & (Freq1 > 0.45 & Freq1 < 0.55))
c <- subset(clumped, index & !(Freq1 > 0.45 & Freq1 < 0.55))
d <- subset(clumped, !index & !(Freq1 > 0.45 & Freq1 < 0.55))

data.frame(
what=c("P C", "NP C", "P NC", "NP NC"),
count=c(nrow(a),
nrow(b),
nrow(c),
nrow(d)),
n=c(mean(a$TotalSampleSize),
mean(b$TotalSampleSize),
mean(c$TotalSampleSize),
mean(d$TotalSampleSize)),
tausq=c(mean(a$tausq),
mean(b$tausq),
mean(c$tausq),
mean(d$tausq))
)

library(data.table)
res <- fread("zcat results/16/16_1.txt.gz")


index <- with(res, 
	((Allele1 == "a" & Allele2 == "t") |
	(Allele1 == "t" & Allele2 == "a") |
	(Allele1 == "g" & Allele2 == "c") |
	(Allele1 == "c" & Allele2 == "g"))
)
table(index)

a <- subset(res, index & (Freq1 > 0.45 & Freq1 < 0.55))
b <- subset(res, !index & (Freq1 > 0.45 & Freq1 < 0.55))
c <- subset(res, index & !(Freq1 > 0.45 & Freq1 < 0.55))
d <- subset(res, !index & !(Freq1 > 0.45 & Freq1 < 0.55))

data.frame(
what=c("P C", "NP C", "P NC", "NP NC"),
count=c(nrow(a),
nrow(b),
nrow(c),
nrow(d)),
n=c(mean(a$TotalSampleSize),
mean(b$TotalSampleSize),
mean(c$TotalSampleSize),
mean(d$TotalSampleSize)),
tausq=c(mean(a$tausq),
mean(b$tausq),
mean(c$tausq),
mean(d$tausq))
)


# We have the reference dataset
# We have the results file 

a <- scan("/panfs/panasas01/sscm/gh13047/sandpit/godmc/relateds/processed_data/genetic_data/data.easyqc.flipped.SNPs.txt", what="character")
b <- subset(target_adjusted, V2 %in% a)


refbim <- fread("zcat resources/genetics/1kg_phase3_eur_allchrs_polymorphic.recoded.nodup.frq.gz")
refbim <- subset(refbim, cptid %in% target$snp)
table(refbim$cptid %in% target_adjusted$V2)
table(target_adjusted$V2 %in% refbim$cptid)

head(target_adjusted[! target_adjusted$V2 %in% refbim$cptid, ])

target_adjusted <- subset(target_adjusted, V2 %in% refbim$cptid)
stopifnot(all(refbim$cptid == target_adjusted$V2))

same <- refbim$a0 == target_adjusted$V5 & refbim$a1 == target_adjusted$V6
wrong_ea <- refbim$a0 == target_adjusted$V6 & refbim$a1 == target_adjusted$V5
flipped <- ! same & ! wrong_ea
table(flipped)
table(wrong_ea)

head(target_adjusted[flipped,])
head(reference[flipped,])

table(paste(target_adjusted$V5[flipped], target_adjusted$V6[flipped]))
table(paste(target_adjusted$V5[wrong_ea], target_adjusted$V6[wrong_ea]))
table(paste(target_adjusted$V5[same], target_adjusted$V6[same]))



