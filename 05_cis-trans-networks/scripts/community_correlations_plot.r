# Analyse the correlation analysis

library(dplyr)
library(ggplot2)

load("../results/corres.rdata")
load("../results/graph.rdata")
dat <- ungroup(dat)
corres <- ungroup(corres)
a <- inner_join(dat, subset(corres, select=c(creg, tcpg, r)), by=c("creg", "tcpg"))
a$rquant <- cut(sqrt(abs(a$r)), breaks=5)
a$rquant2 <- cut(a$r, breaks=5)

p1 <- ggplot(subset(a, !is.na(r)), aes(x=r, y=waldratio)) + 
geom_point(alpha=0.5, aes(colour=rquant)) + 
stat_smooth(method="lm") +
scale_colour_brewer(type="qual") +
labs(x="Community probe correlations in ALSPAC", y="Wald ratio in GoDMC") +
theme(legend.position="none")
ggsave(plot=p1, file="../images/wald_vs_correlation.pdf")

ggplot(subset(a, !is.na(r)), aes(x=r, y=waldratio)) + 
geom_point(alpha=0.5) + 
stat_smooth(method="lm") +
scale_colour_brewer(type="qual") +
labs(x="Community probe correlations in ALSPAC", y="Wald ratio in GoDMC") +
facet_grid(. ~ rquant2, scale="free")


with(subset(a, abs(r) > 0.0), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.1), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.2), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.3), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.4), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.5), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, abs(r) > 0.6), table(sign(waldratio) == sign(r))/sum(!is.na(r)))

with(subset(a, r > 0.0), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.1), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.2), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.3), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.4), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.5), table(sign(waldratio) == sign(r))/sum(!is.na(r)))
with(subset(a, r > 0.6), table(sign(waldratio) == sign(r))/sum(!is.na(r)))

with(subset(a, r > 0.0), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.1), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.2), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.3), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.4), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.5), table(sign(waldratio) == sign(r)))
with(subset(a, r > 0.6), table(sign(waldratio) == sign(r)))

temp <- data_frame(cor=c(mean(corres$r^2, na.rm=TRUE), mean_cor_perms), label=c("True", rep("Permutation", length(mean_cor_perms))))
temp$rn <- rank(temp$cor)
p2 <- ggplot(temp, aes(x=rn, y=cor)) +
geom_point(aes(colour=label)) +
labs(x="", y="Average correlation between CpGs")
ggsave(p2, file="../images/permuted_correlations.pdf")

