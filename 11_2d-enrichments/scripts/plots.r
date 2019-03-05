library(ggplot2)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(bipartite)
library(rvest)

load("../results/difres/difres0.rdata")
load("../data/annotations.rdata")
load("../data/blood.rdata")

anno$code <- paste0(anno$antibody3, " : ", anno$cellType, ", ", anno$tissue, " : ", anno$treatment )
blood$code <- paste0(blood$antibody3, " : ", blood$cellType, ", ", blood$tissue, " : ", blood$treatment )
encode_tfbs <- subset(anno, !duplicated(code) & collection == "encode_tfbs")$index
blood_tfbs <- subset(blood, !duplicated(code) & collection == "encode_tfbs")$index

temp <- subset(anno, select=c(index, code))[encode_tfbs,] %>% filter(!duplicated(code))
temp2 <- subset(blood, index %in% blood_tfbs, select=c(index, code)) %>% filter(!duplicated(code))

stab <- subset(difres, Var1 %in% encode_tfbs & Var2 %in% encode_tfbs & sddif >= 10)

stab <- merge(stab, temp, by.x="Var1", by.y="index", all.x=TRUE)
names(stab)[names(stab) == "code"] <- "snp_tfbs"
stab <- merge(stab, temp, by.x="Var2", by.y="index", all.x=TRUE)
names(stab)[names(stab) == "code"] <- "cpg_tfbs"

stab <- subset(stab, select=c(snp_tfbs, cpg_tfbs, val, Min., Max., sddif)) %>% arrange(desc(sddif))
write.csv(stab, file="../results/2d_enrichment.csv")

l <- data_frame(
  perm = 0:100,
  max=NA,
  sum20=NA,
  maxtfbs=NA,
  sum20tfbs=NA,
  maxblood=NA,
  sum20blood=NA
)

for(i in 0:100)
{
  message(i)
  load(paste0("../results/difres/difres", i, ".rdata"))
  a <- subset(difres, Var1 %in% encode_tfbs & Var2 %in% encode_tfbs)
  b <- subset(difres, Var1 %in% blood_tfbs & Var2 %in% blood_tfbs)
  l$max[i+1] <- max(difres$sddif)
  l$sum20[i+1] <- sum(difres$sddif>20)
  l$maxtfbs[i+1] <- max(a$sddif)
  l$sum20tfbs[i+1] <- sum(a$sddif>20)
  l$maxblood[i+1] <- max(b$sddif)
  l$sum20blood[i+1] <- sum(b$sddif>20)
}

# What is the 95% FDR for tfbs and blood

load("../results/difres/difres0.rdata")
difres$snpanno <- anno$antibody2[difres$Var1]
difres$cpganno <- anno$antibody2[difres$Var2]
difres$sddif2 <- difres$sddif
difres$sddif2[difres$val < difres$Mean] <- difres$sddif2[difres$val < difres$Mean] * -1

summary(difres$sddif2)
# hist(difres$sddif2)
# All values that are 'depleted' are not significant
a <- subset(difres, Var1 %in% encode_tfbs & Var2 %in% encode_tfbs)
b <- subset(difres, Var1 %in% blood_tfbs & Var2 %in% blood_tfbs)
thresh_tfbs <- l$maxtfbs[-1] %>% sort(decreasing = TRUE) %>% .[5]
thresh_blood <- l$maxblood[-1] %>% sort(decreasing = TRUE) %>% .[5]
sum(a$sddif > thresh_tfbs, na.rm=T)
sum(b$sddif > thresh_blood, na.rm=T)

sum(a$sddif > thresh_tfbs) / length(encode_tfbs)^2
sum(b$sddif > thresh_blood) / length(blood_tfbs)^2

sum(a$sddif2 > thresh_tfbs) / length(encode_tfbs)^2
sum(b$sddif2 > thresh_blood) / length(blood_tfbs)^2

sum(a$sddif > max(l$maxtfbs[-1])) / length(encode_tfbs)^2
sum(b$sddif > max(l$maxblood[-1])) / length(blood_tfbs)^2

arrange(b, sddif2) %>% head


##############




temp <- subset(difres, Var1 %in% blood$index & Var2 %in% blood$index & sddif > 22, select=c(snpanno, cpganno, sddif))
temp <- subset(temp, !is.na(snpanno) & !is.na(cpganno))
temp <- cSplit(temp, "snpanno", sep = ";", direction = "long")
temp <- cSplit(temp, "cpganno", sep = ";", direction = "long")
temp <- group_by(temp, cpganno, snpanno) %>%
 	summarise(sddif = mean(sddif), count=n())


blood2 <- subset(blood, collection == "encode_tfbs")
temp <- subset(difres, Var1 %in% blood2$index & Var2 %in% blood2$index & sddif > 12, select=c(snpanno, cpganno, sddif)) %>% arrange(desc(sddif)) %>% subset(!duplicated(paste(snpanno, cpganno)))
temp <- subset(temp, !is.na(snpanno) & !is.na(cpganno))
temp <- cSplit(temp, "snpanno", sep = ";", direction = "long")
temp <- cSplit(temp, "cpganno", sep = ";", direction = "long")
temp <- group_by(temp, cpganno, snpanno) %>%
 	summarise(sddif = mean(sddif), count=n())
dim(temp)

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

bloodi <- subset(anno, index %in% blood2$index)$index
notbloodi <- subset(anno, index %in% tfbs$index & !index %in% blood2$index)$index

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

sig_cpg_anno <- unique(b$Var2)
sig_snp_anno <- unique(b$Var1)
sig_snp 

load("../data/trans_clumped.rdata")
temp <- subset(clumped, snp %in% subset(snpres, anno %in% sig_snp_anno)$snp & cpg %in% subset(cpgres, anno %in% sig_cpg_anno)$cpg, select=c(snp, cpg))
temp <- merge(temp, snpres[snpres$anno %in% sig_snp_anno,], by="snp")
temp <- merge(temp, cpgres[cpgres$anno %in% sig_cpg_anno,], by="cpg")

temp$code <- paste(temp$anno.x, temp$anno.y)
target_code <- paste(b$Var1, b$Var2)
temp2 <- subset(temp, code %in% target_code)

enriched_mqtls <- subset(temp2, !duplicated(paste(cpg, snp)))


clumped$enriched <- paste(clumped$cpg, clumped$snp) %in% paste(enriched_mqtls$cpg, enriched_mqtls$snp)
save(enriched_mqtls, file="../data/enriched_mqtls.rdata")

summary(glm(as.numeric(enriched) ~ nproxies, clumped, family="binomial"))
summary(glm(as.numeric(enriched) ~ Effect, clumped, family="binomial"))
summary(glm(as.numeric(enriched) ~ Effect + nproxies, clumped, family="binomial"))
summary(glm(as.numeric(enriched) ~ TotalSampleSize + Effect + nproxies, clumped, family="binomial"))
with(clumped, tapply(TotalSampleSize, enriched, mean))

# clumped <- as.data.frame(clumped, stringsAsFactors=FALSE)
# clumped <- merge(clumped, snpres, by="snp")
# clumped <- merge(clumped, cpgres, by="cpg")

temp <- subset(snpres, anno %in% blood2$index)




### Clustering the bipartite graph

library(ggnetwork)
library(igraph)
library(tnet)
library(intergraph)

dw2 <- dw
dw2[dw2 != 0] <- 1
rownames(dw2) <- paste("SNP", rownames(dw2))
colnames(dw2) <- paste("CpG", colnames(dw2))
n <- ggnetwork(dw2, layout = "kamadakawai")
n$vn <- as.character(n$vertex.names)
n$vn <- gsub("SNP ", "", n$vn)
n$vn <- gsub("CpG ", "", n$vn)
ggplot(n, aes(x, y, xend = xend, yend = yend)) + 
  geom_edges(color = "grey50") + 
  geom_nodelabel(data = n[ !grepl("SNP", n$vertex.names), ],
                 aes(label = vn),
                 color = "red", fontface="bold") +
  geom_nodelabel(data = n[ grepl("SNP", n$vertex.names), ],
                 aes(label = vn),
                 color = "steelblue", fontface = "bold") +
  theme_blank()
ggsave("../images/bipartite_network.pdf", width=12, height=12)
ggsave("../images/bipartite_network.png", width=12, height=12)



## examples

# irf1 - ezh2 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501188/ (ezh2 mediates silencing of irf1)
# irf1 - smc3 - cohesin
# irf1 - bcl3 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4133603/ (co-downregulated during inflammation in macrophages)
# irf1 - atf3 - atf3 is a negative regulator of cytokines which induce irf1 https://www.ncbi.nlm.nih.gov/pubmed/7688150 https://www.ncbi.nlm.nih.gov/pubmed/16688168
# irf1 - max
# irf1 - tr4





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



## pull out top hits

load("../results/difres/difres0.rdata")
top <- subset(difres, sddif > 10 & Var1 %in% blood2$index & Var2 %in% blood2$index)
top <- merge(top, subset(anno, select=c(index, antibody2)), by.x="Var1", by.y="index")
top <- merge(top, subset(anno, select=c(index, antibody2)), by.x="Var2", by.y="index")
top <- arrange(top, desc(sddif))