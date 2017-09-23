library(dplyr)
library(igraph)

load("../../results/16/16_clumped.rdata")
load("../results/coloc.rdata")
load("../results/creg_tcpg.rdata")


resdat <- inner_join(dat, res, by=c("creg"="exposure", "tcpg"="outcome")) %>%
	filter(H4 > 0.8)
dim(resdat)

gr <- graph_from_data_frame(resdat, directed=TRUE)
# plot(gr)
wc <- cluster_walktrap(gr, steps=20)
class(wc)
length(wc)

mem <- membership(wc)
save(resdat, gr, wc, mem, file="../results/graph.rdata")

q()

cpg <- rbind(
	data_frame(cpg=resdat$creg, chr=resdat$creg_chr, pos=resdat$creg_pos),
	data_frame(cpg=resdat$tcpg, chr=resdat$tcpg_chr, pos=resdat$tcpg_pos)
) %>% filter(!duplicated(cpg))

snp <- subset(clumped, snp %in% resdat$snp) %>% subset(!duplicated(snp)) %>% select(snp, snpchr, snppos, Freq1, snptype)
snp <- subset(snp, !duplicated(snp), select=-c(cpg, cis))
names(snp) <- c("snp", "chr", "pos", "maf", "type")
snp$maf[snp$maf > 0.5] <- 1 - snp$maf[snp$maf > 0.5]

snp_cpg <- rbind(
	data_frame(snp=tcpg$snp, cpg=tcpg$tcpg, pval=tcpg$tcpg_pval, b=tcpg$tcpg_b),
	data_frame(snp=creg$snp, cpg=creg$creg, pval=creg$creg_pval, b=creg$creg_b)
) %>% filter(snp %in% snp$snp, cpg %in% cpg$cpg)

cpg_cpg <- data_frame(cpg1=resdat$creg, cpg2=resdat$tcpg, waldratio=resdat$waldratio, sign=resdat$sign, method="Wald ratio")

cpg$id <- 1:nrow(cpg)
snp$id <- 1:nrow(snp)
cpg_cpg$cpg1 <- cpg$id[match(cpg_cpg$cpg1, cpg$cpg)]
cpg_cpg$cpg2 <- cpg$id[match(cpg_cpg$cpg2, cpg$cpg)]

snp_cpg$snp <- cpg$id[match(snp_cpg$snp, snp$snp)]
snp_cpg$cpg <- cpg$id[match(snp_cpg$cpg, cpg$cpg)]


names(cpg) <- c("name", "chr", "pos:INT", "cpgId:ID(cpg)")
write.csv(cpg, file="../data/cpg.csv", row.names=FALSE)

names(snp) <- c("name", "chr", "pos:INT", "maf:FLOAT", "type", "snpId:ID(snp)")
write.csv(snp, file="../data/snp.csv", row.names=FALSE)

names(snp_cpg) <- c(":START_ID(snp)", ":END_ID(cpg)", "pval:FLOAT", "b:float")
write.csv(snp_cpg, file="../data/snp_cpg.csv", row.names=FALSE)

names(cpg_cpg) <- c(":START_ID(cpg)", ":END_ID(cpg)", "b", "sign", "method")
write.csv(cpg_cpg, file="../data/cpg_cpg.csv", row.names=FALSE)



