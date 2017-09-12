
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

gsea_wrapper <- function(gene.list, universe, MSigDB)
{
	# Make list of gene sets
	a <- list()
	for(i in 1:length(MSigDB))
	{
		a[[i]] <- data.frame(category = names(MSigDB)[i], gene.set = names(MSigDB[[i]]), set = sapply(MSigDB[[i]], length), stringsAsFactors = FALSE)
	}
	gene.sets <- bind_rows(a)
	for(i in 1:nrow(gene.sets))
	{
		gene.sets$overlap[i] <- sum(gene.list %in% MSigDB[[gene.sets$category[i]]][[gene.sets$gene.set[i]]])
	}
	gene.sets$pval <- phyper(gene.sets$overlap, gene.sets$set, length(universe)-gene.sets$set, length(gene.list), lower.tail=FALSE)
	gene.sets$fdr <- p.adjust(gene.sets$pval, "fdr")
	return(gene.sets)
}


permute_wc <- function(wc)
{
	wc$cpg <- sample(wc$cpg, replace=FALSE)
	return(wc)
}

# Run main analysis
# Save full result
# Run permutations
# Save number of sets exceeding FDR 0.05

suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(igraph, quietly = TRUE))
suppressPackageStartupMessages(library(mygene, quietly = TRUE))


args <- commandArgs(T)
jid <- as.numeric(args[1])

out1 <- paste0("../results/annot", jid, ".rdata")
out2 <- paste0("../results/annot_perm", jid, ".rdata")
if(all(file.exists(c(out1, out2)))) q()

load("../data/entrez_genes.rdata")
load("../results/graph.rdata")
universe <- unique(anno$values)

wc <- membership(wc)
wc <- data.frame(cpg=names(wc), membership=as.numeric(wc), stringsAsFactors = FALSE)

count <- table(wc$membership)
keep <- as.numeric(names(count)[count >= 10])
wc <- subset(wc, membership %in% keep)

memlist <- unique(sort(wc$membership))
mem <- memlist[jid]

message(jid, " : ", mem)

gene.list <- subset(wc, membership==mem)$cpg %>% get_gene_from_cpg
res <- gsea_wrapper(gene.list, universe, MSigDB2)
res$mem <- mem
save(res, file=out1)

message("Running permutations")

nperm <- 1000
perms <- data.frame(
	mem = mem,
	perm = 0:nperm,
	minpval = NA,
	sig = NA,
	fdr = NA
)

perms$minpval[1] <- min(res$pval, na.rm=TRUE)
perms$sig[1] <- sum(res$pval < 0.05, na.rm=TRUE)
perms$fdr[1] <- sum(res$fdr < 0.05, na.rm=TRUE)

for(i in 1:nperm)
{
	message(i)
	gene.list2 <- subset(permute_wc(wc), membership==mem)$cpg %>% get_gene_from_cpg
	b <- gsea_wrapper(gene.list2, universe, MSigDB2)
	perms$minpval[i+1] <- min(b$pval, na.rm=TRUE)
	perms$sig[i+1] <- sum(b$pval < 0.05, na.rm=TRUE)
	perms$fdr[i+1] <- sum(b$fdr < 0.05, na.rm=TRUE)
}

save(perms, file=out2)
