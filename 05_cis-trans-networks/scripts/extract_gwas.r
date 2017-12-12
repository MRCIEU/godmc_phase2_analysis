# For each SNP in the clusters find if they are associated with a particular trait

# Take every cluster with more than 10 CpGs

library(dplyr)
library(data.table)

load("../results/communities.rdata")

comms <- group_by(communities, cluster) %>% summarise(nsnp=length(unique(snp)), ncpg=length(unique(cpg))) %>% arrange(desc(nsnp)) %>% filter(nsnp > 5)

selcom <- subset(communities, cluster %in% comms$cluster) %>% arrange(cluster)

# convert to rsids


snp_1kg <- fread("../../10_mr-cpg-gwas/data/eur.bim.orig")

snp_1kg$c1 <- nchar(snp_1kg$V5)
snp_1kg$c2 <- nchar(snp_1kg$V6)

snp_1kg <- subset(snp_1kg, c1 == 1 & c2 == 1)
snp_1kg$snp <- paste0("chr", snp_1kg$V1, ":", snp_1kg$V4, ":SNP")
snp_1kg <- subset(snp_1kg, !duplicated(snp))

selcom <- merge(selcom, subset(snp_1kg, select=c(V2, snp)), by="snp")

write.table(unique(selcom$V2), "../data/cluster_rsids.txt", row=F, col=F, qu=F)


# fn <- list.files("../../10_mr-cpg-gwas/data/extracted/")
# fn <- fn[grepl("filtered_gwas_\\d", fn)]

fn <- list.files("../../data/gwas", full.names=TRUE) %>% grep("txt.gz$", ., value=TRUE)

l <- list()
for(i in 1:length(fn))
{
	message(i, " ", fn[i])
	cmd <- paste0("zfgrep -wf ../data/cluster_rsids.txt ", fn[i], " > ", fn[i], ".cluster")
	l[[i]] <- read.table(paste0(fn[i], ".cluster"), he=FALSE)
	l$i <- i
	l$fn <- fn[i]
}

save(l, file="../data/cluster_extract.rdata")



