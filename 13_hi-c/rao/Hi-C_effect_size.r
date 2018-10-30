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

# Join data, gen overlap id, gen abs effect size
real2 <- as.data.frame(real2)
counts <- full_join(cpgtrans, real2, "code")
counts$overlap <- ifelse(is.na(counts$interaction), "no", "yes")
table(counts$overlap)
counts$effect_abs <- abs(counts$Effect)

# Plot data
#1. Plot density plot of effect size of the 637 overlaps vs. all other interchrom mQTLs
dat <- counts %>%
  group_by(overlap) %>%
  summarise(effect_median = median(effect_abs), effect_mean = mean(effect_abs), effect_max = max(effect_abs))

x2 <- (dat[2,4]+dat[1,4])/2
x3 <- x2[1,1]

pdf("Rao_effect_overlaps_vs_no_overlaps.pdf")
ggplot(counts, aes(x = effect_abs, color = overlap)) +
  geom_density() + 
  geom_vline(aes(xintercept = dat[1,4]), color = "#00BFC4") +
  geom_vline(aes(xintercept = dat[2,4]), color = "#F8766D") +
  labs(title= "Density of Absolute Effect Size", x = "Absolute Effect Size", y = "Density", subtitle = "mQTLs that overlap Hi-C interactions vs. no overaps") +
  scale_color_discrete(name ="mQTL Overlap Status", labels=c("Non-overlapping mQTLs", "Overlapping mQTLs")) +
  theme(legend.position = "bottom") +
  annotate(geom = "text", x = x3, y = 5, label="Max abs effect") +
  annotate("segment", x = x3, xend = dat$effect_max, y = 5-0.2, yend = 4.5, colour = "black", size = 0.5, arrow = arrow(type = "closed", length = unit(0.20,"cm")))
    
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


# TRUE = cis, FALSE = trans, ambivalent = both cis and trans

pdf("Rao_effect_overlaps_vs_cis_cat.pdf")
ggplot(dat1, aes(x = effect_abs, color = cpg_cis)) +
  geom_density() + 
  geom_vline(aes(xintercept = dat2[1,4]), color = "#00BFC4") +
  geom_vline(aes(xintercept = dat2[2,4]), color = "#F8766D") +
  labs(title= "Density of Absolute Effect Size", x = "Absolute Effect Size", y = "Density", subtitle = "mQTL Hi-C overlaps by cis/trans category") +
  scale_color_discrete(name ="mQTL CpG cis/trans Status", labels=c("CpG trans only", "CpG cis and trans")) +
  theme(legend.position = "bottom") +
  annotate(geom = "text", x = x3, y = 5, label="Max abs effect") +
  annotate("segment", x = x3, xend = dat2$effect_max, y = 5-0.2, yend = 4.5, colour = "black", size = 0.5, arrow = arrow(type = "closed", length = unit(0.20,"cm")))

dev.off()

#F8766D
#00BFC4
#scale_color_manual(values=c("#CC6666", "#9999CC"))
  sink()



