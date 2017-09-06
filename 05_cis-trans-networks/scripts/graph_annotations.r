
get_gene_from_cpg <- function(cpglist, out="symbol", annot=anno)
{
	xli <- subset(anno, ind %in% cpglist)$values %>% unique
	if(out=="id")
	{
		require(mygene)
		xli <- queryMany(xli, scopes="symbol", fields="entrezgene", species="human")
	}
	return(xli)
}

gsea_enrichment <- function(gene.list, gene.set, universe)
{
	gene.list2 <- gene.list[gene.list %in% universe]
	gene.set2 <- gene.set[gene.set %in% universe]
	if(length(gene.set2) > 0)
	{
		noverlap <- sum(gene.list2 %in% gene.set2)
		p <- phyper(noverlap, length(gene.set2), length(universe) - length(gene.set2), length(gene.list2), lower.tail=FALSE)
	} else {
		noverlap <- 0
		p <- NA
	}
	d <- data.frame(
		gene_orig = length(gene.list),
		set_orig = length(gene.set),
		universe = length(universe),
		gene = length(gene.list2),
		set = length(gene.set2),
		overlap = noverlap,
		pval = p
	)
	return(d)
}

gsea_wrapper <- function(gene.list, universe)
{
	require(MSigDB)

	# Make list of gene sets
	nom1 <- names(MSigDB)
	l <- list()
	for(i in 1:length(nom1))
	{
		l[[i]] <- data.frame(category = nom1[i], gene.set = names(MSigDB[[nom1[i]]]), stringsAsFactors = FALSE)
	}
	gene.sets <- bind_rows(l)

	# For each gene set perform enrichment
	out <- gene.sets %>%
		group_by(category, gene.set) %>%
		do({
			x <- .
			gsea_enrichment(gene.list, MSigDB[[x$category]][[x$gene.set]], universe)
		})
	out$fdr <- p.adjust(out$pval, "fdr")
	return(out)
}

## Run



suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(igraph, quietly = TRUE))
suppressPackageStartupMessages(library(mygene, quietly = TRUE))
suppressPackageStartupMessages(library(MSigDB, quietly = TRUE))


args <- commandArgs(T)
jid <- as.numeric(args[1])
num <- as.numeric(args[2])

out <- paste0("../results/annot", jid, ".rdata")
if(file.exists(out)) q()

load("../data/entrez_genes.rdata")
load("../results/graph.rdata")
universe <- unique(masterlist$ind)

first <- (jid - 1) * num + 1
last <- min(jid * num, length(wc))
message("Running ", first, " to ", last)

l <- list()
for(i in first:last)
{
	message("Community ", i, " of ", length(wc))
	gene.list <- wc[[i]] %>% get_gene_from_cpg
	l[[i]] <- gsea_wrapper(gene.list, universe)
	l[[i]]$community <- i
}

res <- bind_rows(l)
save(res, file=out)
