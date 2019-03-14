# Script for Enrichment of Rao Hi-C Overlaps
library(data.table)
library(tidyverse)
library(ggplot2)

sink("enrichments_results_effect_size.txt")

# Set WD
setwd("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load data
real <- loadRData("merged_overlaps/nodups.data.Rdata") # real data overlaps
str(real)
real_count <- nrow(real)

# No Duplicate Codes
real2 <- subset(real, !duplicated(real$code))
str(real2)
real_count2 <- nrow(real2) #637

# Load in CpG cis categories
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped_cpgciscategories_251018.Rdata")
cpgtrans <- subset(data, select = c(snp, cpg, cpg_cis, cis, Effect, pval, snpchr, snppos, cpgchr, cpgpos))
cpgtrans$code <- paste(cpgtrans$cpg, cpgtrans$snp)
# Subset for Pval
cpgtrans <- subset(cpgtrans, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE )) 
cpgtrans <- subset(cpgtrans, cpgtrans$cpgchr != cpgtrans$snpchr) #


# Join data, gen overlap id, gen abs effect size
real2 <- as.data.frame(real2)
counts <- full_join(cpgtrans, real2, "code")
counts$overlap <- ifelse(is.na(counts$interaction), "no", "yes")
table(counts$overlap)
counts$effect_abs <- abs(counts$Effect)

# ttest between groups
print("t-test of abs effect size between overlapping and non-ovelapping interchrom trans mQTLs")
with(counts, t.test(effect_abs~overlap))

#Welch Two Sample t-test
#
#data:  effect_abs by overlap
#t = -1.8223, df = 674.13, p-value = 0.06886
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0254129384  0.0009480281
#sample estimates:
#  mean in group no mean in group yes
#0.2242266         0.2364591


# Plot data
#1. Plot density plot of effect size of the 637 overlaps vs. all other interchrom mQTLs
dat <- counts %>%
  group_by(overlap) %>%
  summarise(effect_median = median(effect_abs), effect_mean = mean(effect_abs), effect_max = max(effect_abs))

x2 <- (dat[2,4]+dat[1,4])/2
x3 <- x2[1,1]

#dat = no(1), yes(2)
pdf("Rao_effect_overlaps_vs_no_overlaps.pdf")
ggplot(counts, aes(x = effect_abs, color = overlap)) +
  geom_density() + 
  #geom_vline(aes(xintercept = dat[1,4]), color = "#F8766D") +
  #geom_vline(aes(xintercept = dat[2,4]), color = "#00BFC4") +
  labs(x = "Absolute mQTL effect Size", y = "Density") +
  scale_color_discrete(name ="mQTL Overlap Status", labels=c("Non-overlapping mQTL", "Overlapping mQTL")) +
  theme(legend.position = "bottom") 
  #annotate(geom = "text", x = x3, y = 5, label="Max abs effect") +
  #annotate("segment", x = x3, xend = dat$effect_max, y = 5-0.2, yend = 4.5, colour = "black", size = 0.5, arrow = arrow(type = "closed", length = unit(0.20,"cm")))
    
dev.off()

#2. Plot density plot of effect size of overlaps that are just trans mQTL CpG vs. those that are both cis and trans
# Data set grouped by overlap and trans vs cis and trans
# Numbers:
dat1 <- subset(counts, overlap=="yes")

dat2 <- dat1 %>%
  group_by(cpg_cis) %>%
  summarise(effect_median = median(effect_abs), effect_mean = mean(effect_abs), effect_max = max(effect_abs))
dat2

x2 <- (dat2[1,4]+dat2[2,4])/2
x3 <- x2[1,1]

#ttest
print("t-test of abs effect size between overlapping trans only mQTLs and trans+cis mQTLs")
with(dat1, t.test(effect_abs~cpg_cis))

#Welch Two Sample t-test
#
#data:  effect_abs by cpg_cis
#t = -4.8181, df = 408.12, p-value = 2.047e-06
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.09645713 -0.04055547
#sample estimates:
#  mean in group ambivalent      mean in group FALSE
#0.2084973                0.2770036

# TRUE = cis, FALSE = trans, ambivalent = both cis and trans

pdf("Rao_effect_overlaps_vs_cis_cat.pdf")
ggplot() +
  geom_density(data = dat1, aes(x = effect_abs, color = cpg_cis)) + 
  labs(x = "Absolute mQTL effect size", y = "Density") +
  scale_color_discrete(name ="DNAm site annotation", labels=c("cis and trans", "trans only")) +
  theme(legend.position = "bottom")

dev.off()

#0=trans only, 1=cis/trans, 2=non-overlaps
#counts$cat[counts$cpg_cis=="FALSE"] <- "trans only" #260
#counts$cat[counts$cpg_cis=="ambivalent"] <- "cis and trans" #377
#counts$cat[counts$overlap=="no"] <- "non-overlapping trans mQTLs" #17947

counts$cat[counts$cpg_cis=="FALSE"] <- "0" #260
counts$cat[counts$cpg_cis=="ambivalent"] <- "1" #377
counts$cat[counts$overlap=="no"] <- "2" #17947

counts$cat2[counts$cpg_cis=="FALSE"] <- "trans-overlaps"
counts$cat2[counts$overlap=="no"] <- "no-overlap"

# ttest between groups
print("t-test of abs effect size between overlapping trans only mQTLs and non-ovelapping interchrom trans mQTLs")
with(counts, t.test(effect_abs~cat2))

#counts$cat2 <- factor(counts$cat)

# Plot to show the density of trans only overlaps, cis/trans overlaps and all interchrom mQTLs that do not overlap
pdf("Rao_effect_overlaps_vs_cis_cat2.pdf")
ggplot() +
  geom_density(data = counts, aes(x = effect_abs, color = cat)) +
  labs(x = "Absolute mQTL effect size", y = "Density") +
  scale_color_discrete(name ="DNAm site annotation", labels=c("trans only", "cis and trans", "non-overlapping trans mQTLs")) +
  theme(legend.position = "bottom")

dev.off()

#F8766D
#00BFC4
#scale_color_manual(values=c("#CC6666", "#9999CC"))
sink()
