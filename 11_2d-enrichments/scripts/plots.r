library(ggplot2)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(bipartite)

load("../results/difres/difres0.rdata")
load("../data/annotations.rdata")

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

l <- list()
for(i in 1:100)
{
	message(i)
	load(paste0("../results/difres/difres", i, ".rdata"))
	l[[i]] <- data.frame(perm=i, sddif=difres$sddif)
}

load("../results/difres/difres0.rdata")
l <- bind_rows(l)
l <- rbind(l, data.frame(perm=0, sddif=difres$sddif))
l$lab <- "Permutation"
l$lab[l$perm==0] <- "Real"
ggplot(l, aes(x=sddif)) +
geom_density(aes(fill=lab), alpha=0.5) +
xlim(c(0,20)) +
labs(x="Standard deviation units from the mean", fill="")
ggsave("../images/real_vs_perm.pdf")

