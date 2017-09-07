if a cpg is in a community is it more likely 

library(dplyr)
library(igraph)
load("../results/graph.rdata") # "gr"     "mem"    "resdat" "wc"
load("../../10_mr-cpg-gwas/results/cpg-trait.rdata")
cpg_gwas <- res
load("../../06_mr-gwas-cpg/results/tophits_followup.rdata") # res
gwas_cpg <- res
rm(res)
load("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/snps_gwas.rdata") # gwas

# Combine cpg_gwas and gwas_cpg

library(TwoSampleMR)
ao <- available_outcomes()
a <- subset(ao,
	priority == 1 &
	category %in% c("Disease", "Risk factor") &
	population %in% c("European", "Mixed") & 
	! id %in% c(981, 1026, 1080
		# , 1089, 1090, 1031, 1085, 27, 1083, 276, 278, 1059, 1060, 991, 1074, 9, 10, 12, 980, 1013, 294, 982, 965, 966, 985, 986, 280, 286, 1025, 118, 811, 281, 283, 832, 1076, 1077, 1078, 1071, 1072, 967, 815, 32, 968, 755, 837
	) &
	! author %in% c("Cousminer", "Huffman JE")
)

cpg_gwas <- with(subset(cpg_gwas, p < 1e-5 & H4 > 0.5), data_frame(cpg=exposure, trait=outcome))
gwas_cpg <- with(subset(gwas_cpg, method == "Weighted median" & pval < 1e-5), data_frame(cpg=outcome, trait=as.character(exposure)))

gwas2 <- subset(gwas, !is.na(trait))
index <- match(gwas_cpg$trait, gwas2$exposure)




gsea_wrapper <- function(cpg.list, universe, mrres)
{
	# Make list of trait sets
	mrres$trait <- as.character(mrres$trait)
	a <- as.data.frame(table(mrres$trait))
	names(a) <- c("trait", "set")

	for(i in 1:nrow(a))
	{
		a$overlap[i] <- sum(cpg.list %in% subset(mrres, trait == a$trait[i])$cpg)
	}
	a$pval <- phyper(a$overlap, a$set, universe-a$set, length(cpg.list), lower.tail=FALSE)
	a$fdr <- p.adjust(a$pval, "fdr")
	return(a)
}

a <- gsea_wrapper(wc[[1]], 300000, cpg_gwas)
b <- gsea_wrapper(wc[[1]], 300000, gwas_cpg)

overlap
number in mem
number in cpg-trait

for cpg->trait the universe is number of instrumented cpgs
for trait->cpg the universe is number of tested cpgs

