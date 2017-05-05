library(tidyverse)
library(ggthemes)
load("../results/16/16_clumped.rdata")


cisthresh <- c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
transthresh <- c(1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14)
xc <- filter(clumped, cis)
xt <- filter(clumped, !cis)

# Count number of mQTLs

ar <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	ar[i] <- sum(xc$pval < cisthresh[i])
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- sum(xt$pval < transthresh[i])
}
mqtl_counts <- bind_rows(
	tibble(thresh=cisthresh, count=ar, cis='Cis'),
	tibble(thresh=transthresh, count=art, cis='Trans')
)

ggplot(mqtl_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 6 cohorts")
ggsave("../images/mqtl_counts_thresholds.pdf", width=7, height=7)


## Count number of CpGs

ar <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	ar[i] <- length(unique(subset(xc, pval < cisthresh[i])$cpg))
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- length(unique(subset(xt, pval < transthresh[i])$cpg))
}
cpg_counts <- bind_rows(
	tibble(thresh=cisthresh, count=ar, cis='Cis'),
	tibble(thresh=transthresh, count=art, cis='Trans')
)
ggplot(cpg_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="CpGs with at least one mQTL from meta analysis of 6 cohorts")
ggsave("../images/cpg_counts_thresholds.pdf", width=7, height=7)


## Number of independent SNPs per cis and trans

sig <- subset(clumped, (cis & pval < 1e-7) | (!cis & pval < 1e-14))
clump_counts <- group_by(sig, cpg, cis) %>%
	summarise(n=n())

ggplot(clump_counts, aes(x=n)) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per CpG")
ggsave("../images/clump_counts_cpg.pdf", width=7, height=7)


## Number of hits per SNP

clump_counts_snp <- group_by(sig, snp, cis) %>%
	summarise(n=n()) %>%
	group_by(n, cis) %>%
	summarise(count=n())

ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per SNP")
ggsave("../images/clump_counts_snp.pdf", width=7, height=7)



## Estimate Rsq

# What proportion of rsq is cis vs trans

clumped$rsq <- 2 * clumped$Effect^2 * clumped$Freq1 * (1 - clumped$Freq1)

rsq <- filter(clumped, pval < 5e-8) %>%
	group_by(cpg, cis) %>%
	summarise(rsq=sum(rsq))

rsqt <- filter(clumped, pval < 5e-8) %>%
	group_by(cpg) %>%
	summarise(rsq=sum(rsq))
rsq

ggplot(rsq, aes(x=rsq)) +
geom_density(aes(fill=cis), alpha=0.3) +
labs(x="Rsq ")
ggsave(file="../images/rsq_density_cistrans.pdf", width=12, height=8)

ggplot(rsqt, aes(x=rsq)) +
geom_density() +
labs(x="Rsq ")



# Plot cis against trans
temp <- spread(rsq, key=cis, value=rsq)
names(temp) <- c("cpg", "cis", "trans")
temp <- subset(temp, !is.na(cis) & !is.na(trans))
temp$tot <- temp$cis + temp$trans
temp$rat <- temp$cis / temp$tot
ggplot(temp, aes(x=cis, y=trans)) +
geom_point(alpha=0.06) +
labs(x="Rsq in cis", y="Rsq in trans")
ggsave("../images/rsqcis_vs_rsqtrans.png", width=7, height=7)

ggplot(temp, aes(x=rat)) +
geom_density() +
labs(x="Rsq_cis / (Rsq_cis + Rsq_trans)")
ggsave("../images/rsqcis_vs_rsqtrans_ratio.png", width=7, height=7)




mean(rsqt$rsq)

median(subset(rsqt, cpg %in% temp$cpg)$rsq)
median(subset(rsqt, ! cpg %in% temp$cpg)$rsq)

group_by(rsq, cis) %>% summarise(rsq=sum(rsq)/450000)

median()

## Position of CpG relative to SNP

temp1 <- subset(clumped, cpgchr == snpchr)
temp1$posdif <- temp1$snppos - temp1$cpgpos 
ggplot(subset(temp1, abs(posdif) < 2000000), aes(x=posdif)) +
geom_density() +
labs(x="Distance of SNP from CpG")
ggsave("../images/snp_cpg_distance.pdf")

## Plot maf

clumped$maf <- clumped$Freq1
clumped$maf[clumped$maf > 0.5] <- 1 - clumped$maf[clumped$maf > 0.5]
ggplot(clumped, aes(x=maf, y=abs(Effect))) +
geom_point(alpha=0.04) +
scale_x_log10() +
scale_y_log10()
ggsave("../images/maf_vs_beta.png", width=7, height=7)


## Directions

dir_count <- group_by(clumped, Direction) %>%
	summarise(count=n()) %>%
	filter(count > 100)


ggplot(dir_count, aes(x=Direction, y=count)) +
geom_bar(stat="identity") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
ggsave("../images/directions_all.pdf", width=8, height=6)

## hterogeneity by direction

temp3 <- subset(clumped, Direction %in% dir_count$Direction)
ggplot(temp3, aes(x=Direction, y=HetISq)) +
geom_boxplot(fill="red") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


temp4 <- subset(clumped, pval < 1e-14 & Direction %in% dir_count$Direction)
ggplot(temp4, aes(x=Direction, y=HetISq)) +
geom_boxplot(fill="red") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Isq for mQTL with p < 1e-14")
ggsave("../images/isq_direction.pdf", width=10, height=6)



# Conditioanl 

load("../results/16/16_conditional.rdata")

## Number of independent SNPs per cis and trans

sig <- subset(conditional, (cis & pval < 1e-7) | (!cis & pval < 1e-14))
clump_counts <- group_by(sig, cpg, cis) %>%
	summarise(n=n())

p1 <- ggplot(clump_counts, aes(x=n)) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from conditional analysis (p < 1e-7; 1e-14)", y="mQTLs per CpG")
ggsave(p1, file="../images/conditional_counts_cpg.pdf", width=7, height=7)


## Number of hits per SNP

clump_counts_snp <- group_by(sig, snp, cis) %>%
	summarise(n=n()) %>%
	group_by(n, cis) %>%
	summarise(count=n())

p2 <- ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from conditional analysis (p < 1e-7; 1e-14)", y="mQTLs per SNP")
ggsave(p2, file="../images/conditional_counts_snp.pdf", width=7, height=7)


