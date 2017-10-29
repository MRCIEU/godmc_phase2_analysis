library(dplyr)

fn <- read.csv("../../data/gwas/00info.csv")
fn$row <- 1:nrow(fn)

## Collate the MR analyses

l <- list()
for(i in 1:nrow(fn))
{
	message(i)
	f <- paste0("../results/mr_ld/ext/temp/mr_ld/out_", i, ".rdata")
	if(file.exists(f))
	{
		load(f)
		if("maxnsnp" %in% names(res))
		{
			res <- subset(res, select=-c(maxnsnp)) %>% as_data_frame
		} else {
			res <- ungroup(res)
		}
		l[[i]] <- subset(res, Pvalue < 1e-6)
	} else {
		message("Missing")
	}
}

res <- bind_rows(l)
res <- inner_join(res, fn, by=c("outcome"="row"))
dim(res)

## Collate the colocalisation analyses

threshold1 <- 0.05 / (300000 * 698)
threshold2 <- 0.05 / (300000)


load("../results/cpg_trait_coloc.rdata")
cpgtrait <- merge(res, fn, by.x="outcome", by.y="id")

load("../../data/misc/cpg_pos.rdata")
zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())
cpgtrait <- inner_join(cpgtrait, cpgpos, by=c("exposure"="cpg"))
cpgtrait <- filter(cpgtrait, exposure %in% zhou)
cpgtrait$chr <- as.numeric(gsub("chr", "", cpgtrait$cpgchr))

res2 <- subset(cpgtrait, H4 > 0.8)
res2$p[res2$p < 1e-100] <- 1e-100

res2 <- subset(res2, p < threshold2)
# res2 <- merge(res2, subset(ao, select=c(id, category)), by.x="outcome", by.y="id")

load("../../06_mr-gwas-cpg/results/tophits_followup.rdata")


ressig1 <- subset(res, method == "Weighted mode" & pval < threshold1)
ressig2 <- subset(res, method == "Weighted mode" & pval < threshold2)

ressig2$nom <- strsplit(ressig2$exposure, split=" ") %>% sapply(function(x) x[1])

ressig2$category <- as.factor(ressig2$exposure)
levels(ressig2$category) <- c("Risk factor", "Disease", "Disease", "Disease", "Risk factor", "Risk factor", "Risk factor", "Disease", "Disease", "Disease", "Risk factor", "Disease", "Disease", "Disease", "Metabolites", "Risk factor", "Risk factor", "Risk factor")

temp <- rbind(
	data_frame(cpg=res2$exposure, trait=res2$outcome, cpgpos=res2$cpgpos, chr=res2$chr, p=res2$p, what="cpg", category=res2$category),
	data_frame(cpg=ressig2$outcome, trait=ressig2$exposure, cpgpos=ressig2$cpgpos, chr=ressig2$chr, p=ressig2$pval, what="trait", category=ressig2$category)
)
temp$pval <- -log10(temp$p)
temp$pval[temp$what=="trait"] <- temp$pval[temp$what=="trait"] * -1


library(ggplot2)
ggplot(subset(res2), aes(x=cpgpos, y=-log10(p))) +
geom_point(aes(colour=H4)) +
geom_hline(yintercept=-log10(1e-10)) +
facet_grid(. ~ chr, scale="free", space="free") +
ylim(c(0, -log10(min(res2$p)))) +
theme(legend.position="bottom") +
theme(axis.text.x=element_blank()) 


##


library(ggplot2)
library(ggthemes)

ggplot(temp, aes(x=cpgpos, y=pval)) +
geom_point(size=0.3, aes(colour=category)) +
facet_grid(. ~ chr, scale="free", space="free") +
scale_colour_brewer(type="qual") +
theme_tufte() +
scale_y_continuous(breaks=seq(-100, 100, 10)) +
labs(x="CpG position", y="", colour="") +
theme(
	legend.position="bottom",
	axis.text.x=element_blank(),
	axis.ticks.x=element_blank(),
	panel.grid=element_blank(),
	panel.background=element_rect(fill="grey96", linetype="blank"),
	panel.border=element_blank(),
)
ggsave(file="../images/bidirectional_manhattan.pdf", width=14, height=6)


