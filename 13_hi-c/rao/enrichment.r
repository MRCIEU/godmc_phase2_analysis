# Script for Enrichment of Rao Hi-C Overlaps
library(data.table)
library(tidyverse)
library(ggplot2)

sink("enrichments_results.txt")

# Set WD
setwd("/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/")

# Load data

#file_list <- list.files()

files <- lapply(1:1000, function(x){paste0("merged_nonoverlaps/nodups_data_perm_", x, ".Rdata")})

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

perms <- lapply(files, loadRData) # list of permuted overlaps 
real <- loadRData("merged_overlaps/nodups.data.Rdata") # real data overlaps
str(real)
real_count <- nrow(real)

# perform counts
x <- lapply(perms, nrow) # list (df) of count values
y <- do.call(c, x) # vector of count values

y <- as.data.frame(y)

print(paste0("There are: ", sum(y >= real_count) ," perumation datasets that have a higher overlap count than real data"))

#pval
print(paste0("pval for permutations = ", sum(y >= real_count)/length(perms)))
print(paste0("The mean overlap count for permuted datasets is: ", mean(y$y)))
print(paste0("The min overlap count for permuted datasets is: ", min(y$y)))
print(paste0("The max overlap count for permuted datasets is: ", max(y$y)))

#plots of counts
pdf("Rao_count_densityplot.pdf")
ggplot(as.data.frame(y), aes(x = y)) + 
  #geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  geom_vline(aes(xintercept = mean(y)), 
             color = "#FF6666", linetype = "dashed", size = 1) +
  labs(title = "Density of Overlap Counts", x = "Counts", y = "Density", subtitle = "Hi-C Interactions in Permuted and Real Data") +
  geom_vline(aes(xintercept = real_count), 
             color = "cornflowerblue", linetype = "solid", size = 1) +
  scale_x_continuous(limits = c(400, 1200)) +
  annotate(geom = "text", x = mean(y$y + 150), y = 0.010, label = "Permuted") +
  annotate("segment", x = mean(y$y), xend = mean(y$y + 70), y = 0.010, yend = 0.010, colour = "black", size = 1, arrow = arrow(type = "closed", length = unit(0.20,"cm"))) +
  annotate(geom = "text", x = real_count - 150, y = 0.007, label = "Real") +
  annotate("segment", x = real_count, xend = real_count - 110, y = 0.007, yend = 0.007, colour = "black", size = 1, arrow = arrow(type = "closed", length = unit(0.20,"cm")))
dev.off()
  
#Mean number of counts
#table of counts for perms
#plot of counts
#pvals for perms

sink()



