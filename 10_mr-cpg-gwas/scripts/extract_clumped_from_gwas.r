library(data.table)

fn <- read.csv("../../data/gwas/00info.csv")


get_dat <- function(fn, i)
{
	path <- paste0("../../data/gwas/", fn$newfile[i])
	cmd <- paste0("zcat ", path, " | awk '{ if($7 < 0.00001 || NR == 1) { print $0 }}' > ../data/extracted/filtered_gwas_", fn$id[i], ".txt")
	system(cmd)
	cmd <- paste0("fgrep -wf ../data/extracted/mqtl_clumped.txt ../data/extracted/filtered_gwas_", fn$id[i], ".txt > ../data/extracted/filtered_gwas_mqtl_", fn$id[i], ".txt")
	system(cmd)
	# cmd <- paste0("zfgrep -wf ../data/extracted/mqtl_conditional.txt ", path, " | gzip -c > ../data/extracted/filtered_gwas_mqtl_conditional_", fn$id[i], ".txt.gz")
	# system(cmd)

}

# load("../../results/16/16_conditional.rdata")
load("../../results/16/16_clumped.rdata")

# names(conditional)[names(conditional) == "SNP"] <- "snp"

# res1 <- subset(conditional, 
# 	((pJ < 1e-10 & cis) |
# 	(pJ < 1e-14 & !cis)) &
# 	grepl("SNP", snp)
# )

res2 <- subset(clumped, 
	((pval < 1e-10 & cis) |
	(pval < 1e-14 & !cis)) &
	grepl("SNP", snp)
)


snp_1kg <- fread("../data/eur.bim.orig")

snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)

snp_1kg <- subset(snp_1kg, c1 == 1 & c2 == 1)
snp_1kg$snp <- paste0("chr", snp_1kg$V1, ":", snp_1kg$V4, ":SNP")
snp_1kg <- subset(snp_1kg, !duplicated(snp))

# res1 <- merge(res1, subset(snp_1kg, select=c(V2, snp)), by="snp")
res2 <- merge(res2, subset(snp_1kg, select=c(V2, snp)), by="snp")
dir.create("../data/extracted/", recursive = TRUE)
# write.table(res1$V2, file="../data/extracted/mqtl_conditional.txt", row=F, col=F, qu=F)
write.table(res2$V2, file="../data/extracted/mqtl_clumped.txt", row=F, col=F, qu=F)

library(doParallel)
(no_cores <- detectCores() - 1)
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores, type="FORK")
result <- parLapply(cl, 1:nrow(fn), function(i)
{
	get_dat(fn, i)
})
stopCluster(cl)

for(i in 1:nrow(fn))
{
	message(i, ": ", fn$id[i])
	cmd <- paste0("awk -v var='", fn$id[i], "' '{ print var, $1, $2, $3, $4, $5, $6, $7, $8}' ", "../data/extracted/filtered_gwas_mqtl_", fn$id[i], ".txt > ../data/extracted/filtered_gwas_mqtl2_", fn$id[i], ".txt")
	system(cmd)
	cmd <- paste0("mv ../data/extracted/filtered_gwas_mqtl2_", fn$id[i], ".txt ../data/extracted/filtered_gwas_mqtl_", fn$id[i], ".txt")
	system(cmd)
}

cmd <- "cat ../data/extracted/filtered_gwas_mqtl_[0-9]* | gzip -c > ../results/filtered_gwas_mqtl.txt.gz"
system(cmd)

