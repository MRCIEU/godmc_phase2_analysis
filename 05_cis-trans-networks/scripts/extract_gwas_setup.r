library(dplyr)
library(data.table)
library(magrittr)


load("../../10_mr-cpg-gwas/data/snp_1kg.rdata")

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

out <- group_by(entities, chr) %>%
do({
	ents <- .
	message(ents$chr[1])
	temp <- subset(snp_1kg, V1 == ents$chr[1])
	l <- list()
	for(i in 1:nrow(ents))
	{
		message(i, " of ", nrow(ents))
		l[[i]] <- temp %>%
			mutate(posd = abs(V4 - ents$pos[i])) %>%
			arrange(posd) %>% head(1) %$%
			data_frame(snp_chr=V1, snp_pos=V4, snp_rsid=V2, snp_a1 = V5, snp_a2 = V6, snp_name=snp, snp_distance=posd, id=ents$id[i])
	}
	bind_rows(l)
})

entities <- inner_join(entities, out, by="id")

save(entities, file="../data/entity_info.rdata")

