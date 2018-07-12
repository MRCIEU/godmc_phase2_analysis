library(ggplot2)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(bipartite)
library(rvest)

load("../results/difres/difres0.rdata")
load("../data/annotations.rdata")
load("../data/blood.rdata")

temp <- subset(difres, Var1 %in% blood$index & Var2 %in% blood$index & sddif > 22, select=c(snpanno, cpganno, sddif))
temp <- subset(temp, !is.na(snpanno) & !is.na(cpganno))
temp <- cSplit(temp, "snpanno", sep = ";", direction = "long")
temp <- cSplit(temp, "cpganno", sep = ";", direction = "long")
temp <- group_by(temp, cpganno, snpanno) %>%
 	summarise(sddif = mean(sddif), count=n())


temp <- subset(difres, sddif > 22, select=c(snpanno, cpganno, sddif))
temp <- subset(temp, !is.na(snpanno) & !is.na(cpganno))
temp <- cSplit(temp, "snpanno", sep = ";", direction = "long")
temp <- cSplit(temp, "cpganno", sep = ";", direction = "long")
temp <- group_by(temp, cpganno, snpanno) %>%
 	summarise(sddif = mean(sddif), count=n())


# temp2 <- subset(temp, !snpanno %in% cpganno)
# a <- gvisSankey(temp2[,c("snpanno", "cpganno", "sddif")])
# plot(a)

# temp2 <- group_by(temp2, snpanno) %>%

temp <- subset(temp, !grepl("H3K", snpanno) & !grepl("H3K", cpganno))


snpl <- data.frame(snpanno=unique(temp$snpanno))
snpl$snpid <- 1:nrow(snpl)

cpgl <- data.frame(cpganno=unique(temp$cpganno))
cpgl$cpgid <- 1:nrow(cpgl)


temp <- merge(temp, cpgl)
temp <- merge(temp, snpl)


dw <- matrix(0, nrow(snpl), nrow(cpgl))
dw[as.matrix(temp[,c("snpid", "cpgid")])] <- temp[,"sddif"]

rownames(dw) <- snpl$snpanno
colnames(dw) <- cpgl$cpganno
# heatmap(dw)

dw2 <- dw
dw2[dw2 == 0] <- NA
ii <- cut(t(dw2), breaks = seq(min(dw2, na.rm=T), max(dw2, na.rm=T), len = 100), 
          include.lowest = TRUE)
col <- colorRampPalette(c("lightblue", "blue"))(99)[ii]
# col <- grey((dw[index] - min(dw[index])) / max(dw[index] - min(dw[index])))


pdf("../images/bipartite1.pdf", width=12, height=6)
plotweb(dw, text.rot=90,
	col.interaction=col, bor.col.interaction=col,
	y.width.low=0.03,
	y.width.high=0.03
)
dev.off()

bloodi <- subset(anno, tissue == "blood")$index
notbloodi <- anno$index[!anno$index %in% blood]

l <- list()
for(i in 0:100)
{
	message(i)
	load(paste0("../results/difres/difres", i, ".rdata"))
	difres$sddif2 <- difres$sddif
	difres$sddif2[difres$val < difres$Mean] <- difres$sddif2[difres$val < difres$Mean] * -1
	difres$type <- "Other - other"
	difres$type[difres$Var1 %in% bloodi & difres$Var2 %in% bloodi] <- "Blood - blood"
	difres$type[difres$Var1 %in% notbloodi & difres$Var2 %in% notbloodi] <- "Blood - other"
	l[[i+1]] <- data.frame(perm=i, type=difres$type, sddif=difres$sddif, sddif2=difres$sddif2)
}

l <- bind_rows(l)
l$lab <- "Permutation"
l$lab[l$perm==0] <- "Real"

save(l, file="../results/difres_summary.rdata")

p <- ggplot(l, aes(x=sddif)) +
geom_density(aes(fill=lab), alpha=0.5) +
xlim(c(0,20)) +
labs(x="Standard deviation units from the mean", fill="")
ggsave(p, file="../images/real_vs_perm.pdf")

p <- ggplot(l, aes(x=sddif)) +
geom_density(aes(fill=lab), alpha=0.5) +
xlim(c(0,20)) +
labs(x="Standard deviation units from the mean", fill="") +
facet_grid(type ~ ., scale="free_y")
ggsave(p, file="../images/real_vs_perm2.pdf", height=14)

p <- ggplot(l, aes(x=sddif)) +
geom_density(aes(fill=lab, colour=type), alpha=0.5) +
xlim(c(0,20)) +
labs(x="Standard deviation units from the mean", fill="")
ggsave(p, file="../images/real_vs_perm3.pdf")

l$type2 <- as.character(l$type)
l$type2[l$type == "Blood - blood"] <- "Blood"
l$type2[l$type != "Blood - blood"] <- "Other"
p <- ggplot(l, aes(x=sddif2)) +
geom_density(aes(fill=lab, colour=type2), alpha=0.7) +
labs(x="Standard deviation units from the mean", fill="", colour="Tissue") +
scale_fill_brewer(type="qual") +
scale_colour_brewer(type="qual", palette = 2) +
theme(legend.position=c(0.9,0.8)) 
ggsave(p, file="../images/real_vs_perm4.pdf")


l %>% group_by(type2, perm == 0) %>% summarise(n=n(), m=mean(sddif), s=sd(sddif), ma=max(sddif), np=sum(sddif>20))




### Counts

load("../data/trans_clumped.rdata")
clumped <- merge(clumped, snpres, by="snp")
clumped <- merge(clumped, cpgres, by="cpg")

temp <- subset(snpres, anno %in% blood$index)




### Clustering the bipartite graph

library(ggnetwork)
library(igraph)
library(tnet)
library(intergraph)

dw2 <- dw
dw2[dw2 != 0] <- 1
rownames(dw2) <- paste("SNP", rownames(dw2))
colnames(dw2) <- paste("CpG", colnames(dw2))
n <- ggnetwork(dw2)

ggplot(n, aes(x, y, xend = xend, yend = yend)) + 
  geom_edges(color = "grey50") + 
  geom_nodelabel(data = n[ !grepl("SNP", n$vertex.names), ],
                 aes(label = vertex.names),
                 color = "grey50", label.size = NA) +
  geom_nodelabel(data = n[ grepl("SNP", n$vertex.names), ],
                 aes(label = vertex.names),
                 color = "steelblue", fontface = "bold") +
  theme_blank()


n <- dw2 %*% t(dw2) %>%
  graph_from_adjacency_matrix(mode = "undirected", diag = FALSE, weighted = TRUE)


V(n)$group <- cluster_louvain(n) %>%
  membership %>%
  as.character

temp <- as_adjacency_matrix(n, attr = "weight", sparse = FALSE) %>%
  degree_w 
V(n)$degree <- temp[,3]


ggplot(ggnetwork(n, layout = "kamadakawai"),
       aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(alpha = weight)) +
  geom_nodelabel(aes(label = vertex.names, size = degree, color = group)) +
  scale_alpha_continuous(guide = FALSE) +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  scale_size_continuous(range = c(3, 6), guide = FALSE) +
  theme_blank()
