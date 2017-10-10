library(dplyr)

nom <- list.files("../results/out", pattern="rdata$")
nom <- paste0("../results/out/", nom[grepl("fol", nom)])

l <- list()
for(i in 1:length(nom))
{
	load(nom[i])
	message(i)
	message("nrow: ", nrow(res))
	message("ncpg: ", length(unique(res$outcome)))
	message("nexp: ", length(unique(res$exposure)))
	if(nrow(res) > 0)
	{
		res$nom <- nom[i]
		res$chunk <- i
		l[[i]] <- subset(res, pval < 1e-5, select=c(exposure, outcome, method, b, se, pval, nsnp))
	}
}

res <- bind_rows(l)
res$exposure <- as.character(res$exposure)
res <- as_data_frame(res)

load("../../data/misc/cpg_pos.rdata")
zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())
res <- inner_join(res, cpgpos, by=c("outcome"="cpg"))
res <- filter(res, outcome %in% zhou)
res$chr <- as.numeric(gsub("chr", "", res$cpgchr))
save(res,file="../results/tophits_followup.rdata")


group_by(res, method) %>%
	summarise(n1=sum(pval < 1e-5), n2=sum(pval < 1e-10))


# Signifanct = 
# IVW, Mode

threshold1 <- 0.05 / (300000 * 698)
threshold2 <- 0.05 / (300000)

ressig1 <- subset(res, method == "Weighted mode" & pval < threshold1)
ressig2 <- subset(res, method == "Weighted mode" & pval < threshold2)

table(ressig2$exposure)

library(ggplot2)

p1 <- ggplot(ressig2, aes(x=cpgpos, y=-log10(pval))) +
geom_point(aes(colour=exposure)) +
facet_grid(. ~ chr, scale="free", space="free") +
ylim(c(0, -log10(min(ressig2$pval))))
ggsave(plot=p1, file="../images/ressig2_manhattan.pdf", width=15, height=8)