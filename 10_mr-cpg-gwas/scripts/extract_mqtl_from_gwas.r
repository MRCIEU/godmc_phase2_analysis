library(TwoSampleMR)

ids <- scan("ids.txt")
ao <- available_outcomes()
fn <- subset(ao, id %in% ids, select=c(path, filename, trait, sample_size, nsnp, id))
fn <- fn[order(fn$sample_size, decreasing=TRUE), ]
fn <- subset(fn, !duplicated(trait))


get_dat <- function(fn, i)
{
	path <- paste0("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/passqc/", fn$filename[i])

	cmd <- paste0("awk '{ if($7 < 0.00001 || NR == 1) { print $0 }}' ", path, " > filtered_gwas_", fn$id[i], ".txt")
	system(cmd)

	cmd <- paste0("fgrep -wf mqtl.txt filtered_gwas_", fn$id[i], ".txt > filtered_gwas_mqtl_", fn$id[i], ".txt")
	system(cmd)
}

load("../results/16/16_clumped.rdata")

res <- subset(clumped, 
	((pval < 1e-10 & cis) |
	(pval < 1e-14 & !cis)) &
	grepl("SNP", snp)
)

snp_1kg <- fread("eur.bim.orig")

snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)

snp_1kg <- subset(snp_1kg, c1 == 1 & c2 == 1)
snp_1kg$snp <- paste0("chr", snp_1kg$V1, ":", snp_1kg$V4, ":SNP")
snp_1kg <- subset(snp_1kg, !duplicated(snp))

res <- merge(res, subset(snp_1kg, select=c(V2, snp)), by="snp")
write.table(res$V2, file="mqtl.txt", row=F, col=F, qu=F)


for(i in 1:nrow(fn))
{
	message(i, ": ", fn$id[i])
	get_dat(fn, i)
}

for(i in 1:nrow(fn))
{
	message(i, ": ", fn$id[i])
	cmd <- paste0("awk -v var='", fn$id[i], "' '{ print var, $1, $2, $3, $4, $5, $6, $7, $8}' ", "filtered_gwas_mqtl_", fn$id[i], ".txt > filtered_gwas_mqtl2_", fn$id[i], ".txt")
	# system(cmd)
	cmd <- paste0("mv filtered_gwas_mqtl2_", fn$id[i], ".txt filtered_gwas_mqtl_", fn$id[i], ".txt")
	system(cmd)
}

cmd <- "cat filtered_gwas_mqtl* | gzip -c > filtered_gwas_mqtl.txt.gz"
system(cmd)

