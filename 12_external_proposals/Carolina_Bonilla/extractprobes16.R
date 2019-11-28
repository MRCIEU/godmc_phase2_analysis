library(tidyverse)
library(ggthemes)
library(meffil)
library(dplyr)
library(gridExtra)
library(stringr)

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
r<-read.table("cpgs.txt",sep="\t",he=F)
w<-which(clumped$cpg%in%r$V1)

subset<-clumped[w,]
subset<-subset[subset$cis==TRUE,]
length(unique(subset$cpg))
dim(subset)
