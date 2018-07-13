library(igraph)
library(dplyr)
library(purrr)

load("../results/gwas_clusters_nochr6.rdata")
info <- read.csv("../../data/gwas/00info.csv")
info <- data.frame(fn=gsub(".txt.gz", "", info$newfile), id=info$id)

labels <- subset(res, !duplicated(i), select=c(i, fn))
labels$fn <- gsub("../../data/gwas/", "", labels$fn)
labels$fn <- gsub(".txt.gz", "", labels$fn)
dat <- merge(dat, labels, by.x="id", by.y="i")
sign<-0.05/nrow(dat)
g<-grep("metabolites__",dat$fn)
gwa<-dat[-g,]
gwa<-gwa[which(gwa$binom4<sign),]
gwa.clust<-unique(gwa$clust)

load("../results/graph.rdata")

memcount <- group_by(mem, cluster) %>%
	summarise(n=n()) %>%
	arrange(desc(n))

dat <- group_by(dat, creg) %>%
	mutate(snp1=snp[1])


extract_subset <- function(gr, dat, cpg_list)
{
	keep <- subset(dat, creg %in% cpg_list | tcpg %in% cpg_list)
	keep1 <- subset(keep, select=c(creg, tcpg))
	keep1$type <- "CpG"
	keep2 <- data_frame(creg=keep$snp1, tcpg=keep$creg, type="cis")
	keep3 <- data_frame(creg=keep$snp1, tcpg=keep$tcpg, type="trans")
	keep <- bind_rows(keep1, keep2, keep3)
	keep$lab <- paste(keep$creg, keep$tcpg)
	dup <- keep$lab[duplicated(keep$lab)]
	subset(keep, lab %in% dup)
	keep <- subset(keep, !duplicated(lab))
	# keep <- subset(keep, type == "CpG")
	c <- graph_from_data_frame(keep)
	V(c)$what <- "CpG"
	V(c)[grepl("chr", V(c)$name)]$what <- "SNP"
	return(list(gr=c, dat=keep))
}


(temp <- farthest.nodes(gr))
temp <- all_simple_paths(gr, temp$vertices[1], temp$vertices[2])
temp <- extract_subset(gr, dat, names(temp[[1]]))
pdf(file="../images/plot_longest1.pdf")
plot(temp$gr, 
	layout=layout.fruchterman.reingold, 
	vertex.color=as.numeric(as.factor(V(temp$gr)$what)), 
	vertex.size=(as.numeric(as.factor(V(temp$gr)$what))-5)^2, 
	vertex.label=NA, 
	edge.arrow.size=0, 
	edge.color=as.numeric(as.factor(E(temp$gr)$type))
)
dev.off()

for(i in 1:20)
{
	message(i)
	temp <- extract_subset(gr, dat, mem$cpg[mem$cluster == memcount$cluster[i]])
	pdf(paste0("../images/cluster", memcount$cluster[i], ".pdf"))
	plot(temp$gr, 
		layout=layout.fruchterman.reingold, 
		vertex.color=as.numeric(as.factor(V(temp$gr)$what)), 
		vertex.size=(as.numeric(as.factor(V(temp$gr)$what))-5)^2, 
		vertex.label=NA, 
		edge.arrow.size=0, 
		edge.color=as.numeric(as.factor(E(temp$gr)$type))
	)
	dev.off()
}


for(i in 1:length(gwa.clust))
{
	message(i)
	temp <- extract_subset(gr, dat, mem$cpg[mem$cluster == gwa.clust[i]])
	pdf(paste0("../images/cluster", gwa.clust[i], ".pdf"))
	plot(temp$gr, 
		layout=layout.fruchterman.reingold, 
		vertex.color=as.numeric(as.factor(V(temp$gr)$what)), 
		vertex.size=(as.numeric(as.factor(V(temp$gr)$what))-5)^2, 
		vertex.label=NA, 
		edge.arrow.size=0, 
		edge.color=as.numeric(as.factor(E(temp$gr)$type))
	)
	dev.off()
}

temp <- extract_subset(gr, dat, mem$cpg[mem$cluster == 41])
pdf(paste0("../images/cluster", 41, ".pdf"))
plot(temp$gr, 
	layout=layout.fruchterman.reingold, 
	vertex.color=as.numeric(as.factor(V(temp$gr)$what)), 
	vertex.size=(as.numeric(as.factor(V(temp$gr)$what))-5)^2, 
	vertex.label=NA, 
	edge.arrow.size=0, 
	edge.color=as.numeric(as.factor(E(temp$gr)$type))
)
dev.off()


## What do these look like on circos plots??

library(ggbio)
library(GenomicRanges)

hg19chr <- GRanges(seqnames=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chr19", "chr21", "chr22", "chr23"),
	IRanges(rep(1,23),
		c(249250621,243199373,198022430,191154276,180915260,171115067,155270560,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,63025520,59128983,48129895, 51304566, 159138663)
	)
)

extract_subset_granges <- function(gr, dat, cpg_list)
{
	keep <- subset(dat, creg %in% cpg_list | tcpg %in% cpg_list)
	keep <- tidyr::separate(keep, snp, c("snpchr", "snppos", "snptype"), sep=":", remove=FALSE)
	g1 <- GRanges(seqnames=keep$creg_chr, IRanges(keep$creg_pos, keep$creg_pos))
	g2 <- GRanges(seqnames=keep$tcpg_chr, IRanges(keep$tcpg_pos, keep$tcpg_pos))
	temp <- subset(keep, !duplicated(snp))
	g3 <- GRanges(seqnames=temp$snpchr, IRanges(as.numeric(temp$snppos), as.numeric(temp$snppos)))
	g3$y <- 0

	g1$topos <- g2
	return(list(links=g1, snp=g3))
}

temp <- farthest.nodes(gr)
temp <- all_simple_paths(gr, temp$vertices[1], temp$vertices[2])
temp <- extract_subset_granges(gr, dat, names(temp[[1]]))

dev.new()
ggbio() + 
	circle(temp$links, geom = "link", linked.to = "topos", alpha=0.1) +
	circle(temp$snp, geom = "point", aes(y=y)) +
	circle(hg19chr, geom = "ideo", fill = "gray70", trackWidth=0.1)

for(i in 1:20)
{
	message(i)
	temp <- extract_subset_granges(gr, dat, mem$cpg[mem$cluster == memcount$cluster[i]])
	pdf(paste0("../images/cluster_circle", memcount$cluster[i], ".pdf"))
	p <- ggbio() + 
		circle(temp$links, geom = "link", linked.to = "topos", alpha=0.1) +
		circle(temp$snp, geom = "point", aes(y=y)) +
		circle(hg19chr, geom = "ideo", fill = "gray70", trackWidth=0.1)
	print(p)
	dev.off()
}

temp <- extract_subset_granges(gr, dat, mem$cpg[mem$cluster == 41])
pdf(paste0("../images/cluster_circle", 41, ".pdf"))
p <- ggbio() + 
	circle(temp$links, geom = "link", linked.to = "topos", alpha=0.1) +
	circle(temp$snp, geom = "point", aes(y=y)) +
	circle(hg19chr, geom = "ideo", fill = "gray70", trackWidth=0.1)
print(p)
dev.off()


hub <- hub.score(gr)$vector
auth <- authority.score(gr)$vector

pdf("auth_hub.pdf")
plot(log(hub), log(auth))
dev.off()

pdf(file="../images/plot_all.pdf")
plot(gr, layout=layout.drl, vertex.color=abs(as.numeric(hub > 1)-2), vertex.frame.color=abs(as.numeric(hub > 1)-2), vertex.size = hub, vertex.label = NA, edge.arrow.size = 0, edge.width=1, edge.color="black")
dev.off()


temp <- extract_subset_granges(gr, dat, V(grc)$name)
pdf(file="../images/plot_all_circle.pdf")
p <- ggbio() + 
	circle(temp$links, geom = "link", linked.to = "topos", alpha=0.1) +
	circle(temp$snp, geom = "point", aes(y=y)) +
	circle(hg19chr, geom = "ideo", fill = "gray70", trackWidth=0.1)
print(p)
dev.off()



## Remove small communities


# gr <- graph_from_data_frame(dat, directed=TRUE)

# wc <- cluster_walktrap(as.undirected(gr), steps=1000)
# fg <- cluster_fast_greedy(as.undirected(gr))

# # Simplifying doesn't work:
# gr2a <- igraph::simplify(gr2)
# is_simple(gr2)
# is_simple(gr2a)

# gr2adf <- igraph::as_data_frame(gr2a)

# wc <- cluster_walktrap(gr, steps=20)
# class(wc)
# length(wc)

# mem <- membership(wc)
# mem <- data_frame(cpg=names(mem), cluster=as.numeric(mem))



com <- as_data_frame(mem)
names(com) <- c("cpg", "mem")
a <- group_by(com, mem) %>%
	summarise(n=n()) %>%
	filter(n > 20)

comkeep <- subset(com, mem %in% a$mem)
remcpg <- com$cpg[!com$cpg %in% comkeep$cpg]
grc <- delete_vertices(gr, remcpg)


# plot(grc)
hub <- (hub.score(grc)$vector + 2)^2 - 1
ind <- match(comkeep$cpg, V(grc)$name)
all(V(grc)$name == comkeep$cpg)
V(grc)$community <- comkeep$mem
V(grc)$community2 <- V(grc)$community
V(grc)$community2[duplicated(V(grc)$community2)] <- NA
pdf(file="../images/plot_nosmall.pdf")
plot(grc, layout=layout.fruchterman.reingold, vertex.color=V(grc)$community, vertex.frame.color=V(grc)$community, vertex.size = 2, vertex.label = V(grc)$community2, edge.arrow.size = 0)
dev.off()


## Only keep large communities

com <- as_data_frame(mem)
names(com) <- c("cpg", "mem")
a <- group_by(com, mem) %>%
	summarise(n=n()) %>%
	arrange(desc(n)) %>% head(5)

comkeep <- subset(com, mem %in% a$mem)
remcpg <- com$cpg[!com$cpg %in% comkeep$cpg]
grc <- delete_vertices(gr, remcpg)


hub <- (hub.score(grc)$vector + 2)^2 - 1
ind <- match(comkeep$cpg, V(grc)$name)
all(V(grc)$name == comkeep$cpg)
V(grc)$community <- comkeep$mem
pdf(file="../images/plot_5_largest.pdf")
plot(grc, layout=layout.fruchterman.reingold, vertex.color=V(grc)$community, vertex.frame.color=V(grc)$community, vertex.size = hub, vertex.label = NA, edge.arrow.size = 0)
dev.off()



# SNPs point to CpGs
# CpGs that are directly linked to the same 


library(pathfindR)

table(mem$cpg %in% anno$ind)

out <- list()
for(j in 2:10)
{
	i <- memcount$cluster[j]
	message(i)
	dir.create(paste0("../results/pathfindR/", i), recursive=TRUE)
	setwd(paste0("../results/pathfindR/", i))

	temp1 <- subset(anno, ind %in% subset(mem, cluster == i)$cpg) %>%
		filter(!duplicated(values), values != "")
	temp1$ind <- 5
	temp1$pval <- 1e-8
	out[[j]] <- run_pathfindR(temp1, n_processes=1)
	out$cluster <- i
	setwd("../../../../../scripts")

}

