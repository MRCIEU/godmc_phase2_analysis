# options(dplyr.show_progress = FALSE)

library(TwoSampleMR)
library(MendelianRandomization)
library(dplyr)
library(data.table)

main <- function()
{
	args <- commandArgs(T)
	jid <- as.numeric(args[1])

	message("jid: ", jid)
	dir.create("../results/mr_ld2", show=FALSE)
	out <- paste0("../results/mr_ld2/out_", jid, ".rdata")
	if(file.exists(out)) q()

	message("Reading gwas")
	fn <- read.csv("../../data/gwas/00info.csv")
	extnom <- paste0("../data/extracted/filtered_gwas_mqtl_conditional_ready_", fn$id[param$gwas], ".rdata")

	nsplit <- 10
	param <- expand.grid(split=1:10, gwas=1:nrow(fn))
	param <- param[jid,]
	print(param)
	load(extnom)

	message("Loading conditional")
	load("../data/conditional.rdata")
	conditional <- subset(conditional, SNP %in% ext$SNP)
	cpgcount <- group_by(conditional, exposure) %>%
		summarise(n=n())
	rem <- cpgcount$exposure[cpgcount$n == 1]
	conditional <- subset(conditional, ! exposure %in% rem)
	cpgcount <- subset(cpgcount, n != 1)

	splits <- split(unique(cpgcount$exposure), 1:nsplit)
	# splits[[nsplit]] <- NULL
	sapply(splits, length)

	conditional <- subset(conditional, exposure %in% splits[[param$split]])

	ext$id.outcome <- jid
	ext$outcome <- fn$trait[jid]

	bfileo <- "../../data/ref/condsnps.txt"
	scratch <- paste0("../scratch/temp", jid)
	bfile <- scratch

	write.table(unique(ext$SNP), file=paste0(scratch, ".snplist"), row=F, col=F, qu=F)
	cmd <- paste0(
	"plink ",
	" --bfile ", bfileo,
	" --extract ", paste0(scratch, ".snplist"),
	" --make-bed ", 
	" --out ", scratch
	)
	system(cmd)

	ncpg <- length(splits[[param$split]])
	res <- group_by(conditional, exposure) %>%
		do({
			message(.$exposure[1], " of ", ncpg)
			analyse_one(., ext, jid, bfile)
		})

	system(paste0("rm ", scratch, "*"))
	save(res, file=out)

}

analyse_one <- function(x, ext, jid, bfile)
{
	# message(x$exposure[1], " : ", nrow(x))
	maxsnp <- nrow(x)
	ext_out <- subset(ext, snp %in% x$snp)

	if(nrow(ext_out) <= 1)
	{
		# message("Not enough extracted")
		return(data.frame())
	}

	# Convert SNPs in ext_out
	# x <- subset(x, snp %in% ext_out$snp)
	# index <- match(ext_out$snp, x$snp)
	# stopifnot(all(ext_out$snp == x$snp[index]))
	# ext_out$SNP <- x$SNP[index]
	dat <- try(suppressMessages(harmonise_data(x, ext_out)))

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

	ivw <- try(suppressMessages(MendelianRandomization::mr_ivw(mri)))
	if(nrow(dat) > 2)
	{
		egger <- try(suppressMessages(MendelianRandomization::mr_egger(mri)))
	} else {
		egger <- NULL
	}

	if(!is.null(egger) & class(ivw) != "try-error" & class(egger) != "try-error")
	{
		res <- data.frame(
			outcome=dat$id.outcome[1],
			method=c("IVW", "Egger"),
			Estimate=c(ivw@Estimate, egger@Estimate),
			StdError=c(ivw@StdError, egger@StdError.Est),
			Pvalue=c(ivw@Pvalue, egger@Pvalue.Est),
			heter=c(ivw@Heter.Stat[2], egger@Heter.Stat[2]),
			nsnp=ivw$SNPs
		)
	} else if(class(ivw) != "try-error") {
		res <- data.frame(
			outcome=dat$id.outcome[1],
			method=c("IVW"),
			Estimate=c(ivw@Estimate),
			StdError=c(ivw@StdError),
			Pvalue=c(ivw@Pvalue),
			heter=c(ivw@Heter.Stat[2]),
			nsnp=ivw$SNPs
		)
	} else {
		return(data.frame())
	}
	return(res)
}


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

main()
