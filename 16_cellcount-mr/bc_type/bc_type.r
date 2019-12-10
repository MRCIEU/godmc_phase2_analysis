library(ieugwasr)
gi <- gwasinfo()
bcm <- subset(gi, grepl("Astle", author))
bcmm <- read.table("bc_type.txt", he=F, sep="\t", stringsAsFactors=FALSE)
names(bcmm) <- c("trait", "type")
table(bcm$trait %in% bcmm$trait)
bcm$trait[! bcm$trait %in% bcmm$trait]
bcm <- dplyr::inner_join(bcm, bcmm, by="trait")
save(bcm, file="bc_type.rdata")
