library(dplyr)
library(data.table)
library(magrittr)


load("../../10_mr-cpg-gwas/data/snp_1kg.rdata")
names(snp_1kg) <- c("snp_chr", "snp_rsid", "snp_gp", "snp_pos", "snp_a1", "snp_a2", "snp_c1", "snp_c2", "snp_name")

load("../../results/16/16_clumped.rdata")

load("../results/graph.rdata")
communities <- merge(dat, mem, by.x="creg", by.y="cpg")
communities$creg_chr <- gsub("chr", "", communities$creg_chr)
communities$tcpg_chr <- gsub("chr", "", communities$tcpg_chr)



entities <- bind_rows(
	data_frame(name = communities$creg, chr=gsub("chr", "", communities$creg_chr), pos=communities$creg_pos, cluster=communities$cluster, type="creg_cpg"),
	data_frame(name = communities$tcpg, chr=gsub("chr", "", communities$tcpg_chr), pos=communities$tcpg_pos, cluster=communities$cluster, type="tcpg_cpg"),
	data_frame(name = communities$snp,
		chr=do.call(c, lapply(communities$snp, function(x) strsplit(x, split=":")[[1]][1] %>% gsub("chr", "", .))),
		pos=do.call(c, lapply(communities$snp, function(x) strsplit(x, split=":")[[1]][[2]] %>% as.numeric)),
		cluster=communities$cluster, type="snp")
) %>%
	filter(!duplicated(paste(name, cluster))) %>%
	arrange(as.numeric(chr))
entities$id <- 1:nrow(entities)

table(
	subset(entities, type == "tcpg_cpg") %$% unique(name) %in%
	subset(clumped, cis)$cpg
)
cis <- subset(clumped, cis) %>% arrange(pval) %>% filter(!duplicated(cpg)) %>% ungroup %>% select(cpg, snp)

mentities <- left_join(entities, cis, by=c("name"="cpg"))
mentities <- left_join(mentities, snp_1kg, c("snp"="snp_name"))
names(mentities)[names(mentities) == "snp"] <- "snp_name"

unmatched <- subset(mentities, is.na(snp_rsid))

out <- group_by(unmatched, chr) %>%
do({
	ents <- .
	message(ents$chr[1])
	temp <- subset(snp_1kg, snp_chr == ents$chr[1])
	l <- list()
	for(i in 1:nrow(ents))
	{
		message(i, " of ", nrow(ents))
		l[[i]] <- temp %>%
			mutate(posd = abs(snp_pos - ents$pos[i])) %>%
			arrange(posd) %>% head(1) %$%
			data_frame(snp_chr=snp_chr, snp_pos=snp_pos, snp_rsid=snp_rsid, snp_a1 = snp_a1, snp_a2 = snp_a2, snp_name=snp_name, snp_distance=posd, id=ents$id[i])
	}
	bind_rows(l)
})
outt <- ungroup(out) %>% select(-c(chr))
eout <- inner_join(entities, outt, by="id")

entities <- bind_rows(eout, subset(mentities, !is.na(snp_rsid)))


load("../data/snpcontrolsets_selection.rdata")

targets <- subset(f.all, SNP %in% entities$snp_name)

MAF, GC_freq, CpG_freq, nproxies, tssdist

load("../../results/16/16_clumped.rdata")
load("../../10_mr-cpg-gwas/data/snp_1kg.rdata")
library(dplyr)
snp_1kg <- data_frame(SNP=snp_1kg$snp, a1=snp_1kg$V5, a2=snp_1kg$V6, snp_rsid = snp_1kg$V2)
clumped <- subset(clumped, !snp %in% entities$snp_name & cis)
background <- subset(f.all, SNP %in% clumped$snp)
targets$what <- "target"
background$what <- "background"
genomicinfo <- rbind(targets, background)
genomicinfo <- merge(genomicinfo, snp_1kg, by="SNP")

save(entities, genomicinfo, file="../data/entity_info.rdata")

