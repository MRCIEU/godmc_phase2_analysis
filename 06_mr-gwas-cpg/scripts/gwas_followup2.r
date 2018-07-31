library(dplyr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(parallel)
library(ggplot2)
library(ggrepel)

extract_from_17 <- function(fn, togrep)
{
	require(data.table)
	require(tidyr)
	message("Extracting")
	tmp <- tempfile() %>% basename
	outlist <- paste0(tmp, ".snplist")
	fnn <- paste0(tmp)
	write.table(togrep, file=outlist, row=F, col=F, qu=F)
	cmd <- paste0("zfgrep -f ", outlist, " ", fn, " > ", fnn)
	system(cmd)
	message("Reading")
	b <- fread(fnn)
	unlink(fnn)
	unlink(outlist)
	nom <- as.character(unlist(read.table(fn, nrows=1)))
	names(b) <- nom
	b$Pvalue <- as.numeric(b$Pvalue)
	b <- separate(b, MarkerName, c("snp", "cpg"), "_")
	return(b)
}

do_mr <- function(a, b, gw, chunk)
{
	exposure <- subset(a, exposure == gw)
	exposure$SNP <- exposure$id
	# exposure$exposure <- gw
	exposure <- subset(exposure, select=-c(id, trait))

	outcome <- subset(b, snp %in% exposure$SNP)
	if(nrow(outcome) > 0)
	{
		message("Formatting")
		outcome <- format_data(outcome, type="outcome",
			phenotype_col="cpg",
			snp_col="snp",
			beta_col="Effect",
			se_col="StdErr",
			effect_allele_col="Allele1",
			other_allele_col="Allele2",
			eaf_col="Freq1",
			pval_col="P-value",
			samplesize_col="TotalSampleSize"
		)
		outcome$id.outcome <- outcome$outcome
		exposure$SNP <- tolower(exposure$SNP)
		message("Harmonising")
		dat <- suppressMessages(harmonise_data(exposure, outcome, action=2))
		message(nrow(dat))
		if(nrow(dat) > 0)
		{
			message("Analysing")
			mres <- suppressMessages(mr(dat, metho=c("mr_ivw", "mr_wald_ratio")))
			mres$chunk <- chunk
			return(list(mres=mres, dat=dat))
		} else {
			return(NULL)
		}
	} else {
		return(NULL)
	}
}


mr_followup <- function(dat, chunk, threshold)
{
	message(length(unique(dat$outcome)))
	if(nrow(dat) > 1)
	{
		out <- mr(dat)
		out$chunk <- chunk
		temp <- subset(out, method == "Weighted mode" & pval < threshold)
		if(nrow(temp) > 0)
		{
			message("Plotting ", nrow(temp), " sig hits")
			code <- paste(temp$id.exposure, temp$id.outcome)
			res <- subset(out, paste(id.exposure, id.outcome) %in% code)
			dat <- subset(dat, paste(id.exposure, id.outcome) %in% code)
			p <- mr_scatter_plot(res, dat)
			for(i in 1:length(p))
			{
				ggsave(
					plot=p[[i]] + geom_label_repel(aes(label=gene.exposure), point.padding = unit(0.5, "lines")),
					file=paste0("../results/sig_dat/images/", names(p)[i], "_lab.pdf")
				)
				ggsave(
					plot=p[[i]],
					file=paste0("../results/sig_dat/images/", names(p)[i], ".pdf")
				)
			}
			save(dat, file=paste0("../results/sig_dat/", dat$id.exposure[1], ".", chunk, ".rdata"))
		}
		return(out)
	} else {
		return(NULL)
	}
}

#############

dir.create("../results/sig_dat/", show=F)
dir.create("../results/sig_dat/images", show=F)
load("../data/snps_gwas.rdata")

param <- expand.grid(
	gwas=unique(a$exposure),
	chunk=1:300
)
param$gwas <- as.character(param$gwas)
args <- commandArgs(T)
jid <- as.numeric(args[1])
num <- as.numeric(args[2])
threshold <- as.numeric(args[3])
first <- (jid - 1) * num + 1
last <- min(jid * num, nrow(param))
param <- param[first:last, ]
out <- paste0("../results/out/out", jid, ".rdata")
outfu <- paste0("../results/out/out_fu", jid, ".rdata")
outtemp <- paste0(out, ".temp")
if(file.exists(out)) 
{
	message("Already exists")
	q()
}
if(file.exists(outtemp))
{
	load(outtemp)
	i1 <- length(m) + 1
	j1 <- length(mf) + 1
} else {
	i1 <- 1
	j1 <- 1
}

##############

curr <- "0"
m <- list()
mf <- list()
for(i in i1:nrow(param))
{
	message(i, " of ", nrow(param))
	fn <- paste0("../../results/17/17_", param$chunk[i], ".txt.gz")
	if(fn != curr)
	{
		b <- extract_from_17(fn, unique(a$id))
		curr <- fn
	} else {
		message("Same cpg file")
	}
	out <- do_mr(a, b, param$gwas[i], param$chunk[i])
	m[[i]] <- out$mres

	if(!is.null(out$mres))
	{
		sig <- subset(out$mres, pval < threshold)
		if(nrow(sig) > 0)
		{
			code <- paste(sig$id.exposure, sig$id.outcome)
			mf[[j1]] <- mr_followup(
				subset(out$dat, paste(id.exposure, id.outcome) %in% code), 
				param$chunk[i],
				0.2
			)
			j1 <- j1 + 1
		} 
	}
	save(m, mf, file=outtemp)
}

res <- bind_rows(m)
res_fu <- bind_rows(mf)
save(res, file=out)
save(res_fu, file=outfu)
unlink(outtemp)
