library(tidyverse)
library(ggthemes)
library(meffil)
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

# Count number of mQTLs

arr <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	arr[i] <- sum(xc$PvalueRandom < cisthresh[i])
}
arrt <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	arrt[i] <- sum(xt$PvalueRandom < transthresh[i])
}
mqtl_countsr <- bind_rows(
	tibble(thresh=cisthresh, count=arr, cis='Cis'),
	tibble(thresh=transthresh, count=arrt, cis='Trans')
)




##
y<-meffil.get.features("450k")
ncpg<-length(unique(y$name))


####
#n_independent_regions <- 1000000 used in Frank Dudbridge paper

n_independent_regions <- 1000000
#3 billion basepairs residing in 23 pairs of chromosomes
n_bases <- 3000000000

#SNP-CpG distance is 1 Mb 
cis_window <- 2000000

n_independent_regions_cis <- n_independent_regions / (n_bases / cis_window)
n_independent_regions_trans <- n_independent_regions - n_independent_regions_cis

#number of analysed CpGs - all probe on 450k array
ncpg <- ncpg

ntest_cis <- ncpg * n_independent_regions_cis
ntest_trans <- ncpg * n_independent_regions_trans

# we used pval threshold of 1e-05
exp_cis <- ntest_cis * 1e-4
exp_trans <- ntest_trans * 1e-14
exp_cis
exp_trans
###

p1 <- ggplot(mqtl_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 13 cohorts")
ggsave(plot=p1, file="../images/mqtl_counts_thresholds.pdf", width=7, height=7)

p1 <- ggplot(mqtl_countsr, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="Clumped mQTLs from meta analysis of 13 cohorts")
ggsave(plot=p1, file="../images/mqtl_counts_thresholds_random.pdf", width=7, height=7)


## Overlaps between fixed and random

arf <- rep(0, length(cisthresh))
arr <- rep(0, length(cisthresh))
arb <- rep(0, length(cisthresh))
for(i in 1:length(cisthresh))
{
	arf[i] <- sum(xc$pval < cisthresh[i] & xc$PvalueRandom > cisthresh[i])
	arr[i] <- sum(xc$pval > cisthresh[i] & xc$PvalueRandom < cisthresh[i])
	arb[i] <- sum(xc$pval < cisthresh[i] & xc$PvalueRandom < cisthresh[i])
}
art <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	art[i] <- sum(xt$pval < transthresh[i])
}

artf <- rep(0, length(transthresh))
artr <- rep(0, length(transthresh))
artb <- rep(0, length(transthresh))
for(i in 1:length(transthresh))
{
	artf[i] <- sum(xt$pval < transthresh[i] & xt$PvalueRandom > transthresh[i])
	artr[i] <- sum(xt$pval > transthresh[i] & xt$PvalueRandom < transthresh[i])
	artb[i] <- sum(xt$pval < transthresh[i] & xt$PvalueRandom < transthresh[i])
}


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
p1 <- ggplot(cpg_counts, aes(x=as.factor(-log10(thresh)), y=count)) +
geom_bar(stat="identity") +
geom_text(aes(label=count, y=count+10000), size=3) +
facet_grid(. ~ cis, space="free_x", scale="free_x") +
theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
labs(x="-log10 p threshold", y="CpGs with at least one mQTL from meta analysis of 13 cohorts")
ggsave(plot=p1, file="../images/cpg_counts_thresholds.pdf", width=7, height=7)


## Number of independent SNPs per cis and trans

sig <- subset(clumped, (cis & pval < 1e-7) | (!cis & pval < 1e-14))
clump_counts <- group_by(sig, cpg, cis) %>%
	dplyr::summarise(n=n())

p1 <- ggplot(clump_counts, aes(x=as.factor(n))) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per CpG")
ggsave(plot=p1, file="../images/clump_counts_cpg.pdf", width=7, height=7)


## Number of hits per SNP

clump_counts_snp <- group_by(sig, snp, cis) %>%
	dplyr::summarise(n=n()) %>%
	group_by(n, cis) %>%
	dplyr::summarise(count=n())

p1 <- ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from clumping (p < 1e-7; 1e-14)", y="mQTLs per SNP")
ggsave(plot=p1, file="../images/clump_counts_snp.pdf", width=10, height=7)



## Estimate Rsq

# What proportion of rsq is cis vs trans

clumped$rsq <- 2 * clumped$Effect^2 * clumped$Freq1 * (1 - clumped$Freq1)

rsq <- filter(clumped, pval < 5e-8) %>%
	group_by(cpg, cis) %>%
	dplyr::summarise(rsq=sum(rsq))

rsqt <- filter(clumped, pval < 5e-8) %>%
	group_by(cpg) %>%
	dplyr::summarise(rsq=sum(rsq))
rsq

p1 <- ggplot(rsq, aes(x=rsq)) +
geom_density(aes(fill=cis), alpha=0.3) +
labs(x="Rsq ")
ggsave(plot=p1, file="../images/rsq_density_cistrans.pdf", width=12, height=8)

# ggplot(rsqt, aes(x=rsq)) +
# geom_density() +
# labs(x="Rsq ")



# Plot cis against trans
# temp <- spread(rsq, key=cis, value=rsq)
# names(temp) <- c("cpg", "cis", "trans")
# temp <- subset(temp, !is.na(cis) & !is.na(trans))
# temp$tot <- temp$cis + temp$trans
# temp$rat <- temp$cis / temp$tot
# ggplot(temp, aes(x=cis, y=trans)) +
# geom_point(alpha=0.06) +
# labs(x="Rsq in cis", y="Rsq in trans")
# ggsave("../images/rsqcis_vs_rsqtrans.png", width=7, height=7)

# ggplot(temp, aes(x=rat)) +
# geom_density() +
# labs(x="Rsq_cis / (Rsq_cis + Rsq_trans)")
# ggsave("../images/rsqcis_vs_rsqtrans_ratio.png", width=7, height=7)




# mean(rsqt$rsq)

# median(subset(rsqt, cpg %in% temp$cpg)$rsq)
# median(subset(rsqt, ! cpg %in% temp$cpg)$rsq)

# group_by(rsq, cis) %>% summarise(rsq=sum(rsq)/450000)

#median()

## Position of CpG relative to SNP

temp1 <- subset(clumped, cpgchr == snpchr)
temp1$posdif <- temp1$snppos - temp1$cpgpos 
cisdist<-temp1[which(temp1$cis=="TRUE"),]
mediandist<-round(median(abs(cisdist$posdif))/1000,0)

p1 <- ggplot(subset(temp1, abs(posdif) < 2000000), aes(x=posdif)) +
	geom_density() +
	labs(x="Distance of SNP from CpG") +
	annotate(geom="text", x=0, y=4e-5, label=paste("Median distance = ",mediandist,"kb",sep=""), color="black")
ggsave(plot=p1, file="../images/snp_cpg_distance.pdf")

## Plot maf

clumped$maf <- clumped$Freq1
clumped$maf[clumped$maf > 0.5] <- 1 - clumped$maf[clumped$maf > 0.5]
p1 <- ggplot(clumped, aes(x=maf, y=abs(Effect))) +
geom_point(alpha=0.04) +
scale_x_log10() +
scale_y_log10()
ggsave(plot=p1, file="../images/maf_vs_beta.png", width=7, height=7)


## Directions

dir_count <- group_by(clumped, Direction) %>%
	dplyr::summarise(count=n()) %>%
	filter(count > 100)


p1 <- ggplot(dir_count, aes(x=Direction, y=count)) +
geom_bar(stat="identity") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
ggsave(plot=p1, file="../images/directions_all.pdf", width=8, height=6)

## hterogeneity by direction

# temp3 <- subset(clumped, Direction %in% dir_count$Direction)
# ggplot(temp3, aes(x=Direction, y=-log10(HetPVal))) +
# geom_boxplot(fill="red") +
# theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

# temp33 <- group_by(temp3, Direction) %>% summarise(p = sum(HetPVal < 0.01)/n())

# ggplot(temp33, aes(x=Direction, y=p)) + geom_point() +
# theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))



temp4 <- subset(clumped, pval < 1e-14 & Direction %in% dir_count$Direction)
p1 <- ggplot(temp4, aes(x=Direction, y=HetISq)) +
geom_boxplot(fill="red") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
labs(y="Isq for mQTL with p < 1e-14")
ggsave(plot=p1, file="../images/isq_direction.pdf", width=10, height=6)



# Conditioanl 

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_conditional.rdata")

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


