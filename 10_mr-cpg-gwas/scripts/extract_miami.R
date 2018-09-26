ao <- TwoSampleMR::available_outcomes()
library(RNeo4j)
library(dplyr)
library(magrittr)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- subset(ann, select=c(chr, pos))
ann$chr[ann$chr == "chrY"] <- "chr24"
ann$chr[ann$chr == "chrX"] <- "chr23"
ann$chr <- as.numeric(gsub("chr", "", ann$chr))
ann$cpg <- rownames(ann)


##



miami_trait <- function(trait, graph, ann)
{
	trait <- "Type 2 diabetes"

	l <- list()
	for(i in 1:length(trait))
	{
		message(trait[i])
		b <- cypher(graph,
			paste0("MATCH (n1:cpg)-[r:mr]->(n2) where n2.name = '", trait[i], "' RETURN n1.name, n2.name, r.pval, r.beta")
		)
		b <- b %>% mutate(what="cpg to trait")

		c <- cypher(graph,
			paste0("MATCH (n1)<-[r:mr]-(n2) where n2.name = '", trait[i], "' RETURN n1.name, n2.name, r.pval, r.beta")
		)
		c <- c %>% mutate(what="trait to cpg")

		dat <- bind_rows(b, c) %>% as_data_frame
		dat <- inner_join(dat, as_data_frame(ann), c("n1.name" = "cpg"))
		dat <- dat %$% data_frame(cpg=n1.name, trait=n2.name, pval2 = -log10(r.pval), what, chr, cpgpos=pos, sig=FALSE, chrcol=chr %% 2 == 0)

		dat$pval2[dat$what == "trait to cpg"] <- dat$pval2[dat$what == "trait to cpg"] * -1
		l[[i]] <- dat
	}

	dat <- bind_rows(l)

	load("../results/organised_tophits.rdata")
	temp2 <- temp2[temp2$trait %in% trait, ] %>% select(pval2, chrcol, chr, cpgpos, cpg, trait, what, sig)

	dat <- bind_rows(dat, temp2)
	dat$pval2[dat$pval2 > -log10(1e-60)] <- -log10(1e-60)

	dat21 <- subset(dat, what == "cpg to trait" & !sig) %>% sample_n(5000)

	ind <- which(dat$chr == 21 & dat$what == "trait to cpg") %>% sample(5000)
	dat21$cpg <- dat$cpg[ind]
	dat21$cpgpos <- dat$cpgpos[ind]
	dat21$chr <- 21
	dat <- bind_rows(dat, dat21)
	dat <- sample_n(dat, nrow(dat)) %>% arrange(chr)

	p1 <- ggplot(dat %>% filter(!sig), aes(x=cpgpos, y=pval2)) +
	geom_point(data=dat %>% filter(chrcol), size=0.2, alpha=0.3, aes(colour=trait)) +
	geom_point(data=dat %>% filter(!chrcol), size=0.2, aes(colour=trait)) +
	geom_point(data=dat %>% filter(sig), size=2, colour="black", alpha=1) +
	geom_point(data=dat %>% filter(sig), size=1, aes(colour=trait), alpha=1) +
	facet_grid(. ~ chr, scale="free", space="free") +
	scale_colour_manual(values=c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd")) +
	scale_y_continuous(breaks=seq(min(dat$pval2, na.rm=TRUE)-5, max(dat$pval2, na.rm=TRUE)+5, 20)) +
	labs(x="CpG position", y="", colour="") +
	geom_hline(yintercept=-log10(0.05/sum(dat$what == "trait to cpg" & dat$trait == dat$trait[1])), linetype="dotted") +
	geom_hline(yintercept=log10(0.05/sum(dat$what == "cpg to trait" & dat$trait == dat$trait[1])), linetype="dotted") +
	geom_hline(yintercept=0, linetype="solid") +
	guides(colour=guide_legend(ncol=5)) +
	theme(
		legend.position="top",
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		panel.grid=element_blank(),
		panel.background=element_rect(fill="white", linetype="blank"),
		panel.spacing=unit(0, "lines"),
		panel.border=element_blank(),
		strip.text=element_blank(),
		strip.background=element_rect(fill="white", linetype="blank")
	)
	return(list(dat=dat, plot=p1))
}

graph <- startGraph("http://shark.epi.bris.ac.uk:27600/db/data", username="neo4j", password="123qwe")
trait <- c("Type 1 diabetes", "Type 2 diabetes", "2hr glucose", "Corrected insulin response", "Fasting glucose", "Fasting insulin", "Fasting proinsulin", "HbA1C", "HOMA-B")

table(trait %in% ao$trait)
ids <- ao[ao$trait %in% trait,]

fn <- list.files("../results/coloc")
m <- list()
for(i in 1:length(fn))
{
	message(i)
	load(paste0("../results/coloc/", fn[i]))
	m[[i]] <- subset(res, outcome %in% ids$id) %>% as_data_frame()
}
m <- bind_rows(m) %>% filter(!is.na(H4), nsnp >= 10, H4 > 0.7)

a <- subset(ids, id %in% m$outcome, select=c(id, trait))
m <- inner_join(m, a, c("outcome"="id"))

diab_coloc <- m
save(diab_coloc, file="../results/diab_coloc.rdata")


out <- miami_trait(trait, graph, ann)
ggsave(out$p1, file="../images/diabetes_for_bas.png", width=20, height=10)
