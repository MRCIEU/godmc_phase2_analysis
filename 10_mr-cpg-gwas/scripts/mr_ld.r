# options(dplyr.show_progress = FALSE)

library(TwoSampleMR)
library(MendelianRandomization)
library(dplyr)
library(data.table)

convert_to_chrpos <- function(snp, snp_1kg=snp_1kg)
{
	index <- match(snp, snp_1kg$V2)
	return(snp_1kg$snp[index])
}

get_ld_matrix <- function(snps, bfile, cpg, jid)
{
	fn <- paste0("../scratch/", cpg, "_", jid)
	write.table(snps, file=paste0(fn, ".snplist"), row=F, col=F, qu=F)

	cmd <- paste0("plink ",
		" --bfile ", bfile,
		" --extract ", paste0(fn, ".snplist"),
		" --r square ",
		" --out ", fn, 
		" 1>/dev/null"
	)
	system(cmd)
	rmat <- as.matrix(read.table(paste0(fn, ".ld")))
	colnames(rmat) <- rownames(rmat) <- snps
	cmd <- paste0("rm ", fn, "*")
	system(cmd)
	return(rmat)
}

fn <- read.csv("../../data/gwas/00info.csv")
load("../data/conditional.rdata")
load("../data/snp_1kg.rdata")
extnom <- paste0("zcat ../data/extracted/filtered_gwas_mqtl_conditional_", fn$id, ".txt.gz")

args <- commandArgs(T)
jid <- as.numeric(args[1])

dir.create("../results/mr_ld", show=FALSE)
out <- paste0("../results/mr_ld/out_", jid, ".rdata")

ext <- fread(extnom[jid])
names(ext) <- c("snp", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "sample_size.outcome")

ext <- subset(ext, snp %in% snp_1kg$V2)
ext$SNP <- convert_to_chrpos(ext$snp, snp_1kg)

ext$id.outcome <- jid
ext$outcome <- fn$trait[jid]

bfileo <- "../../data/ref/condsnps.txt"
bfile <- out

write.table(unique(ext$SNP), file=paste0(out, ".snplist"), row=F, col=F, qu=F)
cmd <- paste0(
"plink ",
" --bfile ", bfileo,
" --extract ", paste0(out, ".snplist"),
" --make-bed ", 
" --out ", out
)
system(cmd)


## TESTING
# x <- subset(conditional, exposure == "cg00052772")

res <- group_by(conditional, exposure) %>%
	do({
		x <- .
		message(x$exposure[1], " : ", nrow(x))
		maxsnp <- nrow(x)
		ext_out <- subset(ext, snp %in% x$snp)

		if(nrow(ext_out) <= 1)
		{
			message("Not enough extracted")
			return(data.frame())
		}

		# Convert SNPs in ext_out
		# x <- subset(x, snp %in% ext_out$snp)
		# index <- match(ext_out$snp, x$snp)
		# stopifnot(all(ext_out$snp == x$snp[index]))
		# ext_out$SNP <- x$SNP[index]
		dat <- try(harmonise_data(x, ext_out))

		if(class(dat) == "try-error")
			return(data.frame())

		if(is.null(dat))
			return(data.frame())

		if(nrow(dat) < 2)
			return(data.frame())

		dat <- tidyr::separate(dat, SNP, c("chr", "pos", "type"), ":", remove=FALSE)
		dat$chr <- as.numeric(gsub("chr", "", dat$chr))
		dat$pos <- as.numeric(dat$pos)
		dat <- arrange(dat, chr, pos)
		rmat <- try(get_ld_matrix(dat$SNP, bfile, dat$exposure[1], jid))
		if(class(rmat) == "try-error")
			return(data.frame())

		mri <- mr_input(bx=dat$beta.exposure, by=dat$beta.outcome, bxse=dat$se.exposure, byse=dat$se.outcome, correlation=rmat, snps=dat$SNP, effect_allele=dat$effect_allele.exposure, other_allele=dat$other_allele.exposure, eaf=dat$eaf.exposure)

		ivw <- try(MendelianRandomization::mr_ivw(mri))
		egger <- try(MendelianRandomization::mr_egger(mri))

		if(!is.null(egger) & class(ivw) != "try-error" & class(egger) != "try-error")
		{
			res <- data.frame(
				method=c("IVW", "Egger"),
				Estimate=c(ivw@Estimate, egger@Estimate),
				StdError=c(ivw@StdError, egger@StdError.Est),
				Pvalue=c(ivw@Pvalue, egger@Pvalue.Est),
				heter=c(ivw@Heter.Stat[2], egger@Heter.Stat[2]),
				nsnp=ivw$SNPs,
				maxnsnp=maxsnp
			)
		} else if(class(ivw) != "try-error") {
			res <- data.frame(
				method=c("IVW"),
				Estimate=c(ivw@Estimate),
				StdError=c(ivw@StdError),
				Pvalue=c(ivw@Pvalue),
				heter=c(ivw@Heter.Stat[2]),
				nsnp=ivw$SNPs,
				maxnsnp=maxsnp
			)
		} else {
			return(data.frame())
		}
		return(res)
	})

system(paste0("rm ", out, ".*"))
save(res, file=out)

