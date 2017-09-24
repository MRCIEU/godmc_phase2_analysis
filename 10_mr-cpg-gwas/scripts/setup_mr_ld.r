library(dplyr)
library(data.table)

load("../../results/16/16_conditional.rdata")
snp_1kg <- fread("../data/eur.bim.orig")
snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)
i1 <- snp_1kg$c1 == 1 & snp_1kg$c2 == 1
snp_1kg$snp[i1] <- paste0("chr", snp_1kg$V1[i1], ":", snp_1kg$V4[i1], ":SNP")
snp_1kg$snp[!i1] <- paste0("chr", snp_1kg$V1[!i1], ":", snp_1kg$V4[!i1], ":INDEL")

head(snp_1kg)

dim(conditional)
conditional <- subset(conditional, SNP %in% snp_1kg$snp)
index <- match(conditional$SNP, snp_1kg$snp)
stopifnot(all(conditional$SNP == snp_1kg$snp[index]))
conditional$snp <- snp_1kg$V2[index]

# Only keep conditional SNPs where the CpG has more than 1 hit
conditional <- group_by(conditional, cpg) %>%
	mutate(nsnp=n())

conditional <- subset(conditional, nsnp > 1 & !grepl("INDEL", SNP) & ((cis & pJ < 1e-10) | (!cis & pJ < 1e-14)))
length(unique(conditional$SNP))

temp <- subset(snp_1kg, select=c(snp, V5, V6))
conditional <- inner_join(conditional, temp, by=c("SNP"="snp"))

i1 <- conditional$refA == conditional$V5
i2 <- conditional$refA == conditional$V6

table(i1)
table(i2)

table(i1 | i2)

conditional <- subset(conditional, (i1 | i2))
conditional <- group_by(conditional, cpg) %>%
	mutate(nsnp=n())
conditional <- subset(conditional, nsnp > 1)

i1 <- conditional$refA == conditional$V5
i2 <- conditional$refA == conditional$V6

conditional$othA <- conditional$V6
conditional$othA[i2] <- conditional$V5[i2]

table(conditional$refA == conditional$othA)

conditional <- conditional %>%
	select(SNP, 
		exposure=cpg,
		effect_allele.exposure=refA, 
		other_allele.exposure=othA,
		eaf.exposure=freq,
		beta.exposure=b,
		se.exposure=se,
		pval.exposure=pJ,
		sample_size.exposure=n,
		snp=snp,
		pvalorig=p,
		cis=cis
	) %>% as.tbl
conditional$id.exposure <- conditional$exposure


length(unique(conditional$SNP))
## Make smaller reference dataset
write.table(unique(conditional$SNP), file="../../data/ref/condsnps.txt", row=F, col=F, qu=F)
cmd <- paste0("plink",
	" --bfile ../../data/ref/eur",
	" --extract ../../data/ref/condsnps.txt",
	" --make-bed",
	" --out ../../data/ref/condsnps.txt"
)
system(cmd)

save(conditional, file="../data/conditional.rdata")
