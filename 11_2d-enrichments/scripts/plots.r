library(ggplot2)
library(dplyr)
library(tidyr)

load("../results/difres/difres0.rdata")
load("../data/annotations.rdata")

anno$antibody[is.na(anno$antibody) & anno$collection == "sheffield_dnase"] <- anno$description[is.na(anno$antibody)]

difres$snpanno <- anno$antibody[difres$Var1]
difres$cpganno <- anno$antibody[difres$Var2]

temp <- subset(difres, sddif > 30, select=c(snpanno, cpganno, sddif))
temp <- subset(temp, !duplicated(paste(cpganno, snpanno)))
dw <- spread(temp, key=cpganno, value=sddif, fill=0)
rownames(dw) <- dw$snpanno
dw <- as.matrix(dw[,-1])
heatmap(dw)

library(bipartite)
plotweb(dw)
pdf("../images/bipartite1.pdf")
plotweb(dw,method="normal", text.rot=90)
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

