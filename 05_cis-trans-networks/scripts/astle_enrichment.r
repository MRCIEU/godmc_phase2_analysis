library(dplyr)
library(data.table)

run <- function(extracted, genomicinfo, entities)
{
	runcoms <- function(i)
	{
		message(i)
		j <- 1
		o2 <- list()
		temp2 <- subset(temp, what == "background" | SNP %in% subset(entities, cluster == clust[i])$snp_rsid)
		if(sum(temp2$w == 1) > 5)
		{
			o2[[j]] <- summary(glm(w ~ test + nproxies + CpG_freq + GC_freq, data= temp2, family="binomial"))$coefficients %>% as.data.frame %>% mutate(background="clumped", cluster=clust[i]) %>% slice(2)
			o2[[j]]$ncase <- sum(temp2$w == 1)
			o2[[j]]$ncontrol <- sum(temp2$w == 0)

			j <- j + 1


			temp3 <- subset(temp, what != "background")
			temp3$w[! temp3$SNP %in% subset(entities, cluster == clust[i])$snp_rsid] <- 0
			if(sum(temp3$w == 1) > 5)
			{
				o2[[j]] <- summary(glm(w ~ test, data= temp3, family="binomial"))$coefficients %>% as.data.frame %>% mutate(background="communities", cluster=clust[i]) %>% slice(2)
				o2[[j]]$ncase <- sum(temp3$w == 1)
				o2[[j]]$ncontrol <- sum(temp3$w == 0)
				j <- j + 1
			}
		}
		return(bind_rows(o2))
	}


	ents <- subset(entities, type == "tcpg_cpg")
	temp <- inner_join(extracted, genomicinfo, by=c("SNP" = "snp_rsid")) %>% subset(what == "background" | SNP %in% ents$snp_rsid)
	temp$w <- 1
	temp$w[temp$what == "background"] <- 0
	temp$test <- -log10(pmax(as.numeric(temp$PVAL), 1e-50))
	o <- summary(glm(w ~ test + nproxies + CpG_freq + GC_freq, data= temp, family="binomial"))$coefficients %>% as.data.frame %>% mutate(background="clumped", cluster=NA) %>% slice(2)
	o$ncase <- sum(temp$w == 1)
	o$ncontrol <- sum(temp$w == 0)
	clust <- entities %>% group_by(cluster) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(cluster)

	o2 <- lapply(clust, runcoms) %>% bind_rows
	return(bind_rows(o, o2))
}

####

args <- commandArgs(T)
astle_name <- args[1]
# /projects/MRC-IEU/research/projects/ieu2/p5/021/working/data/Astle2016/27863252-GCST004629-EFO_0004833-Build37.f.tsv.gz
output <- args[2]


load("../../05_cis-trans-networks/data/entity_info.rdata")

r <- fread(paste0("zcat ", astle_name, " | cut -f 2,11"), header=TRUE)
extracted <- subset(r, variant_id %in% genomicinfo$snp_rsid)
names(extracted) <- c("SNP", "PVAL")

res <- run(extracted, genomicinfo, entities)
save(res, file=output)
