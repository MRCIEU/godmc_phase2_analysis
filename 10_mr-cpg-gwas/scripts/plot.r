library(dplyr)
library(ggplot2)
library(ggthemes)

threshold1 <- 0.05 / (300000 * 698)
threshold2 <- 0.05 / (300000)

fn <- read.csv("../../data/gwas/00info.csv")

# load("../results/mr_ld_tophits.rdata")
load("../results/cpg_trait_coloc.rdata")
load("../../data/misc/cpg_pos.rdata")
# zhou <- scan("../../../godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what=character())
cpgtrait <- merge(res, fn, by.x="outcome", by.y="id")
cpgtrait <- inner_join(cpgtrait, cpgpos, by=c("exposure"="cpg"))
# cpgtrait <- filter(cpgtrait, exposure %in% zhou)
cpgtrait$chr <- as.numeric(gsub("chr", "", cpgtrait$cpgchr))
# res2 <- subset(cpgtrait, H4 > 0.8)
res2 <- cpgtrait
res2$p[res2$p < 1e-100] <- 1e-100
res2 <- subset(res2, p < threshold2)

load("../../results/16/16_clumped.rdata")
load("../results/mr_concordance.rdata")
temp8$code <- paste(temp8$outcome.x, temp8$exposure)
res2$code <- paste(res2$trait, res2$exposure)
res2$conc <- res2$code %in% subset(temp8, same_sign & sig2a)$code


##

load("../../06_mr-gwas-cpg/results/mrbase_sig_mhc_sign.rdata")

mult <- subset(res, nsnp >= 6 & what2 == "all")
mult <- group_by(mult, code, exposure) %>%
	summarise(same_sign=sign(b[1]) == sign(b[2]), sig1=pval[2] < 0.05/nrow(mult)*2, sig2 = pval[2] < 0.05, sig3 = pval[1] < threshold2)
table(mult$same_sign, mult$sig2, mult$sig3)
group_by(mult, exposure) %>% summarise(n=n(), nsig1=sum(sig1), nsig2=sum(sig2), nsig3=sum(sig2 & sig3)) %>% arrange( nsig2) %>% as.data.frame
sum(mult$sig2 & mult$sig3)

ressig1 <- subset(res, (method == "Wald ratio" | method == "Inverse variance weighted") & what2=="all" & pval < threshold2)
ressig1$sig <- ressig1$code %in% subset(mult, sig2 & sig3)$code
ressig1 <- inner_join(ressig1, cpgpos, by=c("outcome"="cpg"))
ressig1$chr <- as.numeric(gsub("chr", "", ressig1$cpgchr))
# ressig2 <- subset(res, method == "Inverse variance weighted" & pval < threshold2)
# ressig2$nom <- strsplit(ressig2$exposure, split=" ") %>% sapply(function(x) x[1])
ressig1$pval[ressig1$pval < 1e-100] <- 1e-100

# ressig2$category <- as.factor(ressig2$exposure)
# levels(ressig2$category) <- c("Risk factor", "Disease", "Disease", "Disease", "Risk factor", "Risk factor", "Risk factor", "Disease", "Disease", "Disease", "Risk factor", "Disease", "Disease", "Disease", "Metabolites", "Risk factor", "Risk factor", "Risk factor")

load("../../06_mr-gwas-cpg/results/tophits_followup.rdata")
res <- subset(res, method %in% c("Weighted mode"))
res$chr <- as.numeric(gsub("chr", "", res$cpgchr))

cpgpos$chr <- as.numeric(gsub("chr", "", cpgpos$cpgchr))
cpgpos <- subset(cpgpos, !is.na(chr) & chr %in% 1:23)
rand1 <- data_frame(cpg=cpgpos$cpg, trait="a", chr=cpgpos$chr, p=10^-runif(nrow(cpgpos), min=0, max=-log10(threshold2)), what="cpg to trait", sig=FALSE, cpgpos=cpgpos$cpgpos)
rand2 <- data_frame(cpg=cpgpos$cpg, trait="b", chr=cpgpos$chr, p=10^-runif(nrow(cpgpos), min=0, max=-log10(threshold2)), what="trait to cpg", sig=FALSE, cpgpos=cpgpos$cpgpos)
temp <- rbind(
	data_frame(cpg=res2$exposure, trait=res2$outcome, cpgpos=res2$cpgpos, chr=res2$chr, p=res2$p, what="cpg to trait", sig=res2$H4 > 0.8 & res2$conc),
	data_frame(cpg=ressig1$outcome, trait=ressig1$exposure, cpgpos=ressig1$cpgpos, chr=ressig1$chr, p=ressig1$pval, what="trait to cpg", sig=ressig1$sig),
	data_frame(
		cpg=res$outcome, trait=res$exposure, cpgpos=res$cpgpos, chr=res$chr, p=res$pval, what="trait to cpg", sig=FALSE
		),
	rand1, rand2
	)

temp$p[temp$p < 1e-60] <- 1e-60
temp$pval <- -log10(temp$p)
temp$chrcol <- ifelse(temp$chr %% 2 == 0, "#333333", "#666666")

# temp$pval[temp$what=="trait"] <- temp$pval[temp$what=="trait"] * -1

ggplot(temp %>% filter(!sig), aes(x=cpgpos, y=pval)) +
geom_point(size=0.2, aes(colour=chrcol)) +
geom_point(data=temp %>% filter(sig), size=2, colour="black", alpha=1) +
geom_point(data=temp %>% filter(sig), size=1, colour="red", alpha=1) +
facet_grid(what ~ chr, scale="free", space="free") +
scale_colour_manual(values=c("#bbbbbb", "#888888")) +
# theme_tufte() +
scale_y_continuous(breaks=seq(0, 100, 20)) +
labs(x="CpG position", y="-log10 p", colour="") +
ylim(0, 60) +
geom_hline(yintercept=-log10(threshold2), linetype="dotted") +
theme(
	legend.position="none",
	axis.text.x=element_blank(),
	axis.ticks.x=element_blank(),
	panel.grid=element_blank(),
	panel.background=element_rect(fill="white", linetype="blank"),
	panel.spacing=unit(0, "lines"),
	panel.border=element_blank(),
	strip.text.x=element_blank(),
	strip.background=element_rect(fill="white", linetype="blank"),
)
ggsave(file="../images/bidirectional_manhattan2.png", width=10, height=6)


