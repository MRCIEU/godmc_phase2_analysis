library(dplyr)

nom <- list.files("../results/out", pattern="rdata$")
nom <- paste0("../results/out/", nom[grepl("fol", nom)])

l <- list()
for(i in 1:length(nom))
{
	load(nom[i])
	message(i)
	message("nrow: ", nrow(res))
	message("ncpg: ", length(unique(res$outcome)))
	message("nexp: ", length(unique(res$exposure)))
	if(nrow(res) > 0)
	{
		res$nom <- nom[i]
		res$chunk <- i
		l[[i]] <- subset(res, pval < 1e-5, select=c(exposure, outcome, method, b, se, pval, nsnp))
	}
}

res <- bind_rows(l)
res$exposure <- as.character(res$exposure)
res <- as_data_frame(res)


cpglist <- list()
for(i in 1:300)
{
	message(i)
	cpglist[[i]] <- data_frame(chunk=i, cpg=scan(paste0("../../../godmc_phase1_analysis/07.snp_cpg_selection/lists/cpglist_", i, ".txt"), what="character"))
}
cpglist <- bind_rows(cpglist)

res <- inner_join(res, cpglist, by=c("outcome"="cpg"))

load("../../data/misc/cpg_pos.rdata")
zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())
res <- inner_join(res, cpgpos, by=c("outcome"="cpg"))
res <- filter(res, outcome %in% zhou)
res$chr <- as.numeric(gsub("chr", "", res$cpgchr))
save(res,file="../results/tophits_followup.rdata")


group_by(res, method) %>%
	summarise(n1=sum(pval < 1e-5), n2=sum(pval < 1e-10))



# Signifanct = 
# IVW, Mode

threshold1 <- 0.05 / (300000 * 698)
threshold2 <- 0.05 / (300000)

ressig1 <- subset(res, method == "Weighted mode" & pval < threshold1)
ressig2 <- subset(res, method == "Weighted mode" & pval < threshold2)

table(ressig2$exposure)

ressig2$nom <- strsplit(ressig2$exposure, split=" ") %>% sapply(function(x) x[1])

library(ggplot2)

p1 <- ggplot(subset(ressig2, method == "Weighted mode"), aes(x=cpgpos, y=-log10(pval))) +
geom_point(aes(colour=nom)) +
facet_grid(. ~ chr, scale="free", space="free") +
ylim(c(0, -log10(min(ressig2$pval)))) +
scale_colour_brewer(type="qual", palette="Set3") +
theme(legend.position="bottom") +
theme(axis.text.x=element_blank()) 
ggsave(plot=p1, file="../images/ressig2_manhattan.pdf", width=15, height=8)





##

extract_from_17 <- function(fn, togrep)
{
	require(data.table)
	require(tidyr)
	message("Extracting")
	snplist <- unique(a$id)
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


jid <- as.numeric(commandArgs(T)[1])

load("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata")
a <- subset(gwas, mr_keep.exposure == TRUE)
r <- subset(ressig2, chunk == jid)
a1 <- subset(a, select=c(exposure, id))
r <- inner_join(r, a1, "exposure")

library(data.table)
library(tidyr)
library(TwoSampleMR)
temp <- extract_from_17("../../results/17/17_1.txt.gz", paste0(r$id, "_", r$outcome))
outcome <- format_data(temp, type="outcome",
	phenotype_col="cpg",
	snp_col="snp",
	beta_col="Effect",
	se_col="StdErr",
	effect_allele_col="Allele1",
	other_allele_col="Allele2",
	eaf_col="Freq1",
	pval_col="Pvalue",
	samplesize_col="TotalSampleSize"
)

exposure <- subset(a, exposure %in% r$exposure)
exposure$SNP <- tolower(exposure$id)

dat <- harmonise_data(exposure, outcome)

p <- mr_scatter_plot(mr(dat), dat)
pdf("test.pdf")
dev.off()

save(dat)

