library(LOLA)
library(dplyr)
library(data.table)
library(reshape2)
library(GenomicRanges)
# Get communities of interest


tfbsdb <- loadRegionDB("../../data/lola/scratch/ns5bc/resources/regions/LOLACore/hg19")
load("../results/core_communities_cpg_tophits.rdata")
load("../results/communities.rdata")
load("../results/lola_snp_communities.rdata")
load("../data/lola/cpg_granges.rdata")
load("../data/lola/snp_granges.rdata")
load("../results/graph.rdata")
load("../data/grinfo.rdata")

# load("../../11_2d-enrichments/data/annotations.rdata")



##

# 1. Get top Community annotations
# 2. Get CpGs that are in those communities that have those annotations
# 3. Get SNPs that influence those CpGs
# 4. Get annotations for those SNPs

userSet <- "userSet"
temp <- rbind(
	data_frame(cpg=dat$creg, snp=dat$snp, cpg_chr=dat$creg_chr, cpg_pos=dat$creg_pos, cis=TRUE),
	data_frame(cpg=dat$tcpg, snp=dat$snp, cpg_chr=dat$tcpg_chr, cpg_pos=dat$tcpg_pos, cis=FALSE)
) %>%
filter(!duplicated(paste(cpg, snp)))

l <- list()
for(i in 1:nrow(core_communities_cpg_tophits))
{
	message(i)
	l[[i]] <- extractEnrichmentOverlaps(as.data.frame(core_communities_cpg_tophits)[i,], community_cpgs_separate, tfbsdb) %>% as.data.frame
	l[[i]] <- subset(temp, 
			paste(cpg_chr, cpg_pos) %in% paste(l[[i]]$seqnames, l[[i]]$start)
		) %>%
	filter(!duplicated(paste(cpg, snp)))
	l[[i]]$userSet <- core_communities_cpg_tophits$userSet[i]
	l[[i]] <- merge(l[[i]], as.data.frame(core_communities_cpg_tophits)[i,])
}

out <- bind_rows(l)
out$antibody[is.na(out$antibody)] <- out$description[is.na(out$antibody)]

ind <- grepl("CTCF", out$antibody)
out$antibody[ind] <- "CTCF"
ind <- grepl("SMC3", out$antibody)
out$antibody[ind] <- "SMC3"
ind <- grepl("Znf143", out$antibody)
out$antibody[ind] <- "Znf143"
out$antibody <- toupper(out$antibody)


out <- out %>%
	arrange(desc(pValueLog), cpg, antibody) %>%
	group_by(userSet, cpg, antibody) %>%
	mutate(count=n()) %>%
	filter(!duplicated(paste(userSet, cpg, antibody)))



### SNP annotations

xsnps <- subset(grinfo, snp %in% unique(out$snp) & !duplicated(snp))
nom <- xsnps$snp
xsnps <- GRanges(seqnames=xsnps$snpchr, IRanges(xsnps$min, xsnps$max))
names(xsnps) <- nom

snpo <- lapply(tfbsdb$regionGRL, function(x){
	xsnps[queryHits(findOverlaps(xsnps, x))]
})

snpo2 <- list()
for(i in 1:length(snpo))
{
	if(length(snpo[[i]]) > 0)
	{
		x <- snpo[[i]][!duplicated(names(snpo[[i]]))]
		snpo2[[i]] <- as.data.frame(x)
		snpo2[[i]]$ind <- i
	}
}
snpo2 <- bind_rows(snpo2)
snpo2 <- cbind(snpo2, temp[snpo2$ind, ])


# temp <- tfbsdb[[2]]
# subset(temp, collection == "encode_tfbs" & cellType == "HFF")

# m <- list()
# for(i in 163:nrow(out))
# {
# 	message(i)
# 	ind <- which(temp$collection == out$collection[i] & temp$cellType == out$cellType[i])
# 	message(length(ind))
# 	n <- list()
# 	for(j in 1:length(ind))
# 	{
# 		temp2 <- tfbsdb$regionGRL[[ind[j]]]
# 		temp3 <- xsnps[names(xsnps) == out$snp[i]]
# 		temp3[queryHits(findOverlaps(temp3, temp2))]
# 		if(length(temp) > 0)
# 		{
# 			ret <- temp[ind[j], ]
# 			ret$snp <- names(temp3)
# 			n[[j]] <- ret
# 		}
# 	}
# 	m[[i]] <- bind_rows(n)
# }


save(snpo2, out, m, file="../results/plotting_community_enrichments.rdata")


#####



load("../results/plotting_community_enrichments.rdata")


snpo2$antibody[is.na(snpo2$antibody)] <- snpo2$description[is.na(snpo2$antibody)]

snpo2$snp <- paste0(snpo2$seqnames, ":", snpo2$start, ":SNP")

table(snpo2$snp, snpo2$antibody)
table(paste(snpo2$snp, snpo2$antibody)) %>% sort

snpanno <- group_by(snpo2, snp, antibody) %>%
	summarise(n=n()) %>%
	arrange(desc(n)) %>%
	filter(between(row_number(), 1, pmin(10, n())))


tabb <- data_frame(snp=unique(snpanno$snp) %>% sort, b=out$snp %>% unique %>% sort)
snpanno <- inner_join(snpanno, tabb)
snpanno$snp <- b
stopifnot(all(snpanno$b %in% out$snp))
stopifnot(all(out$snp %in% snpanno$b))



bigdat <- list(

	select(snpanno, start=antibody, end=snp, weight=n) %>%
		mutate(section=1),
	filter(out, !duplicated(paste(snp, cpg))) %>%
		ungroup() %>%
		mutate(weight=as.numeric(cis) * -1 + 2, section=2) %>%
		select(start=snp, end=cpg, weight=weight, section=section),
	out %>% ungroup %>%
		mutate(section=3) %>%
		select(start=cpg, end=antibody, weight=pValueLog, section=section)
) %>% bind_rows()

library(googleVis)
bigdat$start[bigdat$section==1] <- paste(bigdat$section[bigdat$section==1], bigdat$start[bigdat$section==1])
bigdat$end[bigdat$section==3] <- paste(bigdat$section[bigdat$section==3], bigdat$end[bigdat$section==3])
a <- gvisSankey(subset(bigdat), from="start", to="end", weight="")
plot(a)


library(networkD3)

nom <- tibble(name=unique(c(bigdat$start, bigdat$end)))
nom$index <- 1:nrow(nom)-1

bigdat <- merge(bigdat, nom, by.x="start", by.y="name")
names(bigdat)[names(bigdat) == "index"] <- "source"
bigdat <- merge(bigdat, nom, by.x="end", by.y="name")
names(bigdat)[names(bigdat) == "index"] <- "target"
bigdat$value <- 1

sankeyNetwork(Links=bigdat, Nodes=nom, Source="source", Target="target", Value="value", NodeID="name")


make_layer <- function(d, label1, label2)
{
	data.frame(input = paste0(label1, "_", d[,1]), output = paste0(label2, "_", as.character(d[,2])), iname = as.character(d[,1]), oname = as.character(d[,2]), stringsAsFactors=FALSE)
}

tabb$snp2 <- tabb$snp
tabb <- subset(tabb, select=-c(snp))
out <- merge(out, tabb, by.x="snp", by.y="b")
out$snp <- out$snp2

temp <- filter(out) %>% 
	ungroup %>%
	mutate(section=3) %>%
	select(start=cpg, end=antibody, weight=pValueLog, section=section) %>% as.data.frame
l3 <- make_layer(temp, "CpG", "Target")
temp <- filter(out, !duplicated(paste(snp, cpg))) %>%
		ungroup() %>%
		mutate(weight=as.numeric(cis) * -1 + 2, section=2) %>%
		select(start=snp, end=cpg, weight=weight, section=section) %>% as.data.frame
l2 <- make_layer(temp, "SNP", "CpG")
temp <- select(snpanno, start=antibody, end=snp, weight=n) %>%
		mutate(section=1) %>% as.data.frame
l1 <- make_layer(temp, "Input", "SNP")




d <- rbind(l1, l2, l3)



g <- graph.data.frame(d)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- name
t <- match(V(g)$type, unique(type) )
# t[t %in% c(1,4)] <- NA
V(g)$name[t %in% 2:3]  <- ""

set.seed(1)
l <- layout_with_fr(g, minx=t, maxx=t, grid="nogrid")
# ind1 <- l[,1] == 4 & l[,2] < 0
# ind2 <- l[,1] == 4 & l[,2] > 0
# l[, 2][ind1] <- l[, 2][ind1] + 3
# l[, 2][ind2] <- l[, 2][ind2] - 11
# # l[, 2][l[,1] == 4] <- l[, 2][l[,1] == 4] * 2
# l[, 2][l[,1] == 4] <- l[, 2][l[,1] == 4] * 8
# l[, 2][l[,1] == 1] <- l[, 2][l[,1] == 1] * 3

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

for(i in unique(l[,1]))
{
	message(i)
	l[,2][l[,1] == i] <- range01(rank(l[,2][l[,1] == i])) * 10
}

# l <- layout_as_bipartite(g, types=V(g)$type == "CpG")
size <- rep(5, length(V(g)$type))
size[V(g)$type == "SNP"] <- 6
size[V(g)$type == "Target"] <- 4
pdf("../images/community_enrichments.pdf")
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=size, vertex.label.degree=-pi/2, label=V(g)$name, vertex.label.cex=0.5, vertex.frame.color=t, vertex.label.color="black")
dev.off()


####


g <- graph.data.frame(d)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- ""
t <- match(V(g)$type, unique(type) )
l <- layout_with_fr(g, miny=t, maxy=t)
l <- layout_as_bipartite(g, types=V(g)$type == "CpG")
size <- rep(1, length(V(g)$type))
size[V(g)$type == "SNP"] <- 2
size[V(g)$type == "Target"] <- 4
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=size)

########


make_layer <- function(d, label1, label2)
{
	data.frame(input = paste0(label1, "_", d[,1]), output = paste0(label2, "_", as.character(d[,2])), iname = as.character(d[,1]), oname = as.character(d[,2]), stringsAsFactors=FALSE)
}


a <- subset(snpres, snp %in% names(community_snps_separate$STAG1) & anno == 667) %>% filter(!duplicated(snp))
b <- subset(communities, snp %in% a$snp)[,1:2]
l1 <- make_layer(as.data.frame(b[,c(2,1)]), "SNP", "CpG")
c <- subset(cpgres, cpg %in% b$cpg & anno %in% c(499, 810, 90, 863))
l2 <- make_layer(as.data.frame(c[, c(2,1)]), "CpG", "Target")
l2 <- subset(l2, !duplicated(paste(input, output)))
l1 <- subset(l1, oname %in% c$cpg)
l0 <- data.frame("CEBPB", subset(a, snp %in% l1$iname)$snp) %>% make_layer("Input", "SNP")
d <- rbind(l0, l1, l2)

g <- graph.data.frame(d)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- name
t <- match(V(g)$type, unique(type) )
# t[t %in% c(1,4)] <- NA
V(g)$name[t %in% 2:3]  <- ""
l <- layout_with_fr(g, minx=t, maxx=t, grid="nogrid")
# l <- layout_as_bipartite(g, types=V(g)$type == "CpG")
size <- rep(5, length(V(g)$type))
size[V(g)$type == "SNP"] <- 6
size[V(g)$type == "Target"] <- 4
pdf("test.pdf")
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=size, vertex.label.degree=-pi/2, label=V(g)$name, vertex.label.cex=0.2, vertex.frame.color=t, vertex.label.color="black")
dev.off()

plot(g, vertex.label.rotate=-pi/2, label=V(g)$name)



g <- graph.data.frame(l0)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- ""

l <- layout_as_bipartite(g, types=V(g)$type == "Input")
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=1)

g <- graph.data.frame(l1)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- ""

l <- layout_as_bipartite(g, types=V(g)$type == "SNP")
png("l1.png")
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=1, edge.color="darkorange")
dev.off()

g <- graph.data.frame(l2)
type <- do.call(rbind, strsplit(V(g)$name, split="_"))[,1]
name <- do.call(rbind, strsplit(V(g)$name, split="_"))[,2]
V(g)$type <- type
V(g)$name <- ""

l <- layout_as_bipartite(g, types=V(g)$type == "CpG")
png("l2.png")
plot(g, vertex.color=t, layout=l, edge.arrow.size=0.2, vertex.size=1, edge.color="chocolate3")

dev.off()

a <- subset(snpres, snp %in% names(community_snps_separate$STAG1) & anno == 667) %>% filter(!duplicated(snp))
b <- subset(communities, snp %in% a$snp)[,1:2]
b$val <- 1
br <- reshape2::dcast(b, snp ~ cpg, fill=0)
rownames(br) <- br$snp
br <- as.matrix(br[,-1])

temp <- as.data.frame(as.table(br), stringsAsFactors=FALSE)


g <- graph.data.frame(b)
V(g)$type <- substr(V(g)$name, 1, 2)
t <- match(V(g)$type, c("cg", "ch") )

l <- layout_with_fr(g, miny=t, maxy=t)
plot(g, vertex.color=t, layout=l, edge.width=E(g)$val)


ar <- subset(ar, snp %in% br$snp)
ar <- reshape2::dcast(ar, snp ~ pert)


d <- data.frame(target=c(90,499,810,863), func="Cohesin")



cr <- reshape2::dcast(c, anno ~ cpg)

1 3 7 9

extract




enr_snp_communities_top <- subset(enr_snp_communities, grepl("encode", collection)) %>% 
	arrange(desc(pValueLog)) %>% 
	mutate(pbonf = p.adjust(10^(-pValueLog), "bonferroni"), pfdr = p.adjust(10^(-pValueLog), "fdr")) %>%
	filter(pfdr < 0.1) %>% 
	data.table
enr_snp_communities_top


snpi <- c(2,4,7,8,9)

l <- list()
for(i in 1:length(snpi))
{
	l[[i]] <- extractEnrichmentOverlaps(enr_snp_communities_top[snpi[i],], community_snps_separate$STAG1, tfbsdb)

}


load("")


```



# Get CpGs in those communities



# Get SNPs that influence those CpGs



g <- sample_pa(20, m=2)
minC <- rep(-Inf, vcount(g))
maxC <- rep(Inf, vcount(g))
minC[1] <- maxC[1] <- 0
co <- layout_with_fr(g, minx=minC, maxx=maxC,
                               miny=minC, maxy=maxC)
co[1,]
plot(g, layout=co, vertex.size=1:20, edge.arrow.size=0.2,
  vertex.label=c("ego", rep("", vcount(g)-1)), rescale=FALSE,
  xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
  vertex.label.color="red")
axis(1)
axis(2)
