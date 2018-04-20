library(dplyr)
library(rvest)
library(tidyr)



t1 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=1&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t2 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=2&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()
t3 <- html("https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?ra=encode/cv.ra&type=Cell+Line&tier=3&bgcolor=FFFEE8") %>% html_nodes("table") %>% .[[1]] %>% html_table()

encode <- bind_rows(t1, t2, t3)

# Index annotations
anno$index <- 1:nrow(anno)

# Everything without a tissue but with a celltype
temp1 <- subset(anno, is.na(tissue) & !is.na(cellType))
encodet <- subset(encode, select=c(cell, Tissue))
temp1$tissue <- encodet$Tissue[
        match(gsub("-", "", tolower(temp1$cell)), gsub("-", "", tolower(encodet$cell)))
]

load("../results/core_communities_cpg_tophits.rdata")
load("../results/core_global_cpg.rdata")
load("../results/ext_communities_cpg_tophits.rdata")
load("../results/ext_global_cpg.rdata")
load("../results/s25_communities_cpg_tophits.rdata")
load("../results/s25_global_cpg.rdata")

core_global_cpg <- subset(core_global_cpg, collection == "encode_tfbs")
core_global_cpg$fdr <- p.adjust(10^-core_global_cpg$pValueLog, "fdr")

table(encode$cell %in% core_global_cpg$cellType)
table(core_communities_cpg_tophits$cellType %in% encode$cell)

index <- match(core_global_cpg$cellType, encode$cell)
core_global_cpg$tissue <- encode$Tissue[index]

index <- match(core_communities_cpg_tophits$cellType, encode$cell)
table(core_communities_cpg_tophits$cellType == encode$cell[index])
core_communities_cpg_tophits$tissue <- encode$Tissue[index]

states <- read.table("../data/roadmap_states.txt", he=T, sep="\t")
tissues <- read.table("../data/jul2013.roadmapData_tissues.txt", he=T, sep="\t")

s25_global_cpg <- separate(s25_global_cpg, filename, c("ct", "v1", "v2", "v3", "state"), sep="_")
s25_global_cpg$state <- gsub("E", "", s25_global_cpg$state)
s25_global_cpg$state <- gsub(".bed", "", s25_global_cpg$state)
index <- match(s25_global_cpg$state, states$STATE)
s25_global_cpg$antibody <- states$MNEMONIC[index]
s25_global_cpg$cellType <- s25_global_cpg$ct
index <- match(s25_global_cpg$cellType, tissues$ID)
s25_global_cpg$tissue <- tissues$Tissue[index]

s25_communities_cpg_tophits <- separate(s25_communities_cpg_tophits, filename, c("ct", "v1", "v2", "v3", "state"), sep="_")
s25_communities_cpg_tophits$state <- gsub("E", "", s25_communities_cpg_tophits$state)
s25_communities_cpg_tophits$state <- gsub(".bed", "", s25_communities_cpg_tophits$state)
index <- match(s25_communities_cpg_tophits$state, states$STATE)
s25_communities_cpg_tophits$antibody <- states$MNEMONIC[index]
s25_communities_cpg_tophits$cellType <- s25_communities_cpg_tophits$ct
index <- match(s25_communities_cpg_tophits$cellType, tissues$ID)
s25_communities_cpg_tophits$tissue <- tissues$Tissue[index]

rm(index, t1, t2, t3, encode, tissues, states)
save.image("../results/lola_organised.rdata")