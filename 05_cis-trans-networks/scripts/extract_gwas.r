if(!require(gwasvcftools))
{
	if(!required(devtools)) install.packages("devtools")
	devtools::install_github("MRCIEU/gwasvcftools")
}
library(gwasvcftools)
library(argparse)

# create parser object
parser <- ArgumentParser()
parser$add_argument('--entities', required=TRUE)
parser$add_argument('--bcf-dir', required=TRUE)
parser$add_argument('--gwas-id', required=TRUE)
parser$add_argument('--out', required=TRUE)
parser$add_argument('--bfile', required=TRUE)
parser$add_argument('--get-proxies', default='yes')
parser$add_argument('--vcf-ref', required=FALSE)
parser$add_argument('--tag-r2', type="double", default=0.6)
parser$add_argument('--tag-kb', type="double", default=5000)
parser$add_argument('--tag-nsnp', type="double", default=5000)
parser$add_argument('--palindrome-freq', type="double", default=0.4)
parser$add_argument('--threads', type="integer", default=1)
parser$add_argument('--no-clean', action="store_true", default=FALSE)

args <- parser$parse_args()
print(args)
tempname <- tempfile(pattern="extract", tmpdir=dirname(args[['out']]))
bcf <- file.path(args[['bcf_dir']], args[['gwas_id']], "harmonised.bcf")
load(args[['entities']])

# snplist <- subset(entities, type != "creg_cpg") %$%
# 	data_frame(V1=snp_chr, V2=snp_pos, V3=snp_a1, V4=snp_a2, V5="b37", V6=snp_rsid) %>%
# 	subset(!duplicated(paste(V1,V2)))

snplist <- genomicinfo %$%
	data_frame(V1=gsub("chr", "", chromosome), V2=position, V3=a1, V4=a2, V5="b37", V6=snp_rsid) %>%
	subset(!duplicated(paste(V1,V2)))

extracted <- gwasvcftools::extract(bcf, snplist, tempname, args[['get_proxies']], args[["bfile"]], args[["vcf_ref"]], threads=args[["threads"]])
extracted$mrbaseid <- args[['gwas_id']]

temp <- subset(extracted, duplicated(ID))$ID
extracted1 <- subset(extracted, !ID %in% temp)
extracted2 <- subset(extracted, ID %in% temp) %>% filter(is.na(proxy.chrom))
extracted <- rbind(extracted1, extracted2)
save(extracted, file=args[["out"]])


run <- function(extracted, genomicinfo, entities)
{
	runcoms <- function(i)
	{
		message(i)
		j <- 1
		o2 <- list()
		temp2 <- subset(temp, what == "background" | SNP %in% subset(entities, cluster == clust[i])$snp_name)
		if(sum(temp2$w == 1) > 5)
		{
			o2[[j]] <- summary(glm(w ~ test + nproxies + CpG_freq + GC_freq, data= temp2, family="binomial"))$coefficients %>% as.data.frame %>% mutate(background="clumped", cluster=clust[i]) %>% slice(2)
			o2[[j]]$ncase <- sum(temp2$w == 1)
			o2[[j]]$ncontrol <- sum(temp2$w == 0)

			j <- j + 1


			temp3 <- subset(temp, what != "background")
			temp3$w[! temp3$SNP %in% subset(entities, cluster == clust[i])$snp_name] <- 0
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
	genomicinfo$chromosome <- gsub("chr", "", as.character(genomicinfo$chromosome)) %>% as.numeric
	temp <- inner_join(extracted, genomicinfo, by=c("CHROM" = "chromosome", "POS"="position")) %>% subset(what == "background" | SNP %in% ents$snp_name)
	temp$w <- 1
	temp$w[temp$what == "background"] <- 0
	temp$test <- -log10(pmax(temp$PVAL, 1e-50))
	o <- summary(glm(w ~ test + nproxies + CpG_freq + GC_freq, data= temp, family="binomial"))$coefficients %>% as.data.frame %>% mutate(background="clumped", cluster=NA) %>% slice(2)
	o$ncase <- sum(temp$w == 1)
	o$ncontrol <- sum(temp$w == 0)
	clust <- entities %>% group_by(cluster) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(cluster)

	o2 <- lapply(clust, runcoms) %>% bind_rows
	return(bind_rows(o, o2))
}

res <- run(extracted, genomicinfo, entities)
save(res, file=paste0(args[["gwas_id"]], "_enr.rdata"))
