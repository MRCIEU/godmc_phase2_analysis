library(ggplot2)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(bipartite)
library(rvest)

load("../results/difres/difres0.rdata")
load("../data/annotations.rdata")

# Get tissue types for encode

t1 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=1&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t2 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=2&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t3 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=3&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()

encode <- bind_rows(t1, t2, t3)

# Index annotations
anno$index <- 1:nrow(anno)

# Everything without a tissue but with a celltype
temp1 <- subset(anno, is.na(tissue) & !is.na(cellType))
encodet <- subset(encode, select=c(cell, Tissue))
temp1$tissue <- encodet$Tissue[
	match(gsub("-", "", tolower(temp1$cell)), gsub("-", "", tolower(encodet$cell)))
]

table(is.na(temp1$tissue), temp1$collection)

# Merge back

temp2 <- subset(anno, !index %in% temp1$index)
anno2 <- rbind(temp1, temp2)
anno2 <- anno2[order(anno2$index),]
stopifnot(all(anno2$filename == anno$filename))
anno <- anno2

anno$antibody[is.na(anno$antibody) & anno$collection == "sheffield_dnase"] <- anno$description[is.na(anno$antibody) & anno$collection == "sheffield_dnase"]

anno$antibody[is.na(anno$antibody) & anno$collection == "encode_segmentation"] <- anno$description[is.na(anno$antibody) & anno$collection == "encode_segmentation"]

anno$antibody[is.na(anno$antibody) & anno$collection == "ucsc_features"] <- anno$description[is.na(anno$antibody) & anno$collection == "ucsc_features"]


anno$antibody[anno$collection == "sheffield_dnase"] <- NA

anno$antibody2 <- sapply(strsplit(as.character(anno$antibody), split="_"), function(x) x[1])
anno$antibody2 <- sapply(strsplit(as.character(anno$antibody2), split=" "), function(x) x[1])
anno$antibody2 <- sapply(strsplit(as.character(anno$antibody2), split="\\("), function(x) x[1])
anno$antibody2 <- gsub("eGFP-", "", anno$antibody2)


difres$snpanno <- anno$antibody2[difres$Var1]
difres$cpganno <- anno$antibody2[difres$Var2]
difres$sddif2 <- difres$sddif
difres$sddif2[difres$val < difres$Mean] <- difres$sddif2[difres$val < difres$Mean] * -1

summary(difres$sddif2)
# hist(difres$sddif2)
# All values that are 'depleted' are not significant

blood <- subset(anno, tissue == "blood")

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

blood <- subset(anno, tissue == "blood")$index
notblood <- anno$index[!anno$index %in% blood]

l <- list()
for(i in 0:100)
{
	message(i)
	load(paste0("../results/difres/difres", i, ".rdata"))
	difres$sddif2 <- difres$sddif
	difres$sddif2[difres$val < difres$Mean] <- difres$sddif2[difres$val < difres$Mean] * -1
	difres$type <- "Other - other"
	difres$type[difres$Var1 %in% blood & difres$Var2 %in% blood] <- "Blood - blood"
	difres$type[difres$Var1 %in% notblood & difres$Var2 %in% notblood] <- "Blood - other"
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