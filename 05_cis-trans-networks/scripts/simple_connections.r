# Construct c-reg t-cpg networks

library(tidyverse)
library(igraph)

load("../../results/16/16_clumped.rdata")
load("../../data/misc/naeem.rdata")
naeem_keep <- subset(naeem, Flag.discard.keep. == "keep")$probe
table(unique(clumped$cpg) %in% naeem_keep) / length(unique(clumped$cpg))
table(naeem_keep %in% naeem$probe) / nrow(naeem)
clumped <- subset(clumped, cpg %in% naeem_keep)


# Get the trans hits

tcpg <- subset(clumped, !cis & pval < 1e-14)
tcpg <- data_frame(snp=tcpg$snp, tcpg_chr=tcpg$cpgchr, tcpg_pos=tcpg$cpgpos, tcpg_pval=tcpg$pval, tcpg_b=tcpg$Effect, tcpg=tcpg$cpg)

creg <- subset(clumped, cis & pval < 1e-14)
creg <- data_frame(snp=creg$snp, creg_chr=creg$cpgchr, creg_pos=creg$cpgpos, creg_pval=creg$pval, creg_b=creg$Effect, creg=creg$cpg)

dat <- inner_join(creg, tcpg, by="snp")
dat <- select(dat, creg, tcpg, snp, creg_chr, tcpg_chr, creg_pos, tcpg_pos, creg_pval, tcpg_pval, creg_b, tcpg_b)
dat$waldratio <- dat$tcpg_b / dat$creg_b
dat$sign <- sign(dat$waldratio)

dat
table(table(dat$creg))

table(dat$creg_b > dat$tcpg_b)
table(sign(dat$creg_b) == sign(dat$tcpg_b))

gr <- graph_from_data_frame(dat, directed=TRUE)
plot(gr)



pathlengths <- data_frame(cpg=unique(c(dat$creg)), max_length=NA, connection=NA)

for(i in 1:nrow(pathlengths))
{
	message(i)
	temp <- all_simple_paths(gr, from=pathlengths$cpg[i])
	pathlengths$max_length[i] <- max(sapply(temp, length))
	pathlengths$connection[i] <- length(temp)
}

table(pathlengths$max_length)
devtools::install_github("nicolewhite/RNeo4j")

library(RNeo4j)
gr <- startGraph("http://localhost:7474/db/data/", username="neo4j", password="123qwe")
clear(gr, FALSE)
addConstraint(gr, "CpG", "cpg")
addConstraint(gr, "SNP", "snp")
for(i in 1:nrow(dat))
{
	message(i, " of ", nrow(dat))
	a <- getOrCreateNode(gr, "CpG", cpg=dat$creg[i], pos=dat$creg_pos[i], chr=dat$creg_chr[i])
	b <- getOrCreateNode(gr, "CpG", cpg=dat$tcpg[i], pos=dat$tcpg_pos[i], chr=dat$tcpg_chr[i])
	createRel(a, "c", b, pval=dat$tcpg_pval[i], b=dat$waldratio[i], sign=dat$sign[i])
}

dats <- subset(dat, !duplicated(paste(snp, creg)))
for(i in 1:nrow(dats))
{
	message(i, " of ", nrow(dats)
	a <- getOrCreateNode(gr, "CpG", cpg=dats$creg[i], pos=dats$creg_pos[i], chr=dats$creg_chr[i])
	s <- getOrCreateNode(gr, "SNP", snp=dats$snp[i])
	createRel(s, "i", a, pval=dats$creg_pval[i], b=dats$creg_b[i])
}


summary(gr)


## Create CSVs for upload

# CpGs
## - cpg
## - chr
## - pos

# SNPs
## - snp
## - chr
## - pos
## - MAF
## - type

# SNP - CpGs
## - snp
## - cpg
## - b
## - se
## - pval
## - n

# CpG-CpG
## - wald ratio



tcpg <- subset(clumped, !cis & pval < 1e-14)
tcpg <- data_frame(snp=tcpg$snp, tcpg_chr=tcpg$cpgchr, tcpg_pos=tcpg$cpgpos, tcpg_pval=tcpg$pval, tcpg_b=tcpg$Effect, tcpg=tcpg$cpg)

creg <- subset(clumped, cis & pval < 1e-14)
creg <- data_frame(snp=creg$snp, creg_chr=creg$cpgchr, creg_pos=creg$cpgpos, creg_pval=creg$pval, creg_b=creg$Effect, creg=creg$cpg)

dat <- inner_join(creg, tcpg, by="snp")
dat <- select(dat, creg, tcpg, snp, creg_chr, tcpg_chr, creg_pos, tcpg_pos, creg_pval, tcpg_pval, creg_b, tcpg_b)
dat$waldratio <- dat$tcpg_b / dat$creg_b
dat$sign <- sign(dat$waldratio)

cpg <- rbind(
	data_frame(cpg=dat$creg, chr=dat$creg_chr, pos=dat$creg_pos),
	data_frame(cpg=dat$tcpg, chr=dat$tcpg_chr, pos=dat$tcpg_pos)
) %>% filter(!duplicated(cpg))

snp <- subset(clumped, snp %in% dat$snp) %>% subset(!duplicated(snp)) %>% select(snp, snpchr, snppos, Freq1, snptype)
snp <- subset(snp, !duplicated(snp), select=-c(cpg, cis))
names(snp) <- c("snp", "chr", "pos", "maf", "type")
snp$maf[snp$maf > 0.5] <- 1 - snp$maf[snp$maf > 0.5]

snp_cpg <- rbind(
	data_frame(snp=tcpg$snp, cpg=tcpg$tcpg, pval=tcpg$tcpg_pval, b=tcpg$tcpg_b),
	data_frame(snp=creg$snp, cpg=creg$creg, pval=creg$creg_pval, b=creg$creg_b)
) %>% filter(snp %in% snp$snp, cpg %in% cpg$cpg)

cpg_cpg <- data_frame(cpg1=dat$creg, cpg2=dat$tcpg, waldratio=dat$waldratio, sign=dat$sign, method="Wald ratio")

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


