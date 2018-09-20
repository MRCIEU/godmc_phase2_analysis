library(TwoSampleMR)
library(tidyverse)

threshold <- 1.4e-7

load("../data/snps_gwas.rdata")
## loads a
gkeep <- filter(a, grepl("mr", data_source.exposure))
gkeeps <- gkeep %>% group_by(id.exposure, trait, exposure) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup %>% filter(!duplicated(trait))
rm(a)


load("../results/mrbase_tophits_full.rdata")
res <- subset(res, id.exposure %in% gkeeps$id.exposure)
het <- subset(het, id.exposure %in% gkeeps$id.exposure)
plei <- subset(plei, id.exposure %in% gkeeps$id.exposure)
res$code <- paste(res$id.exposure, res$id.outcome)
het$code <- paste(het$id.exposure, het$id.outcome)
het <- subset(het, method == "Inverse variance weighted")
plei$code <- paste(plei$id.exposure, plei$id.outcome)

codes <- data_frame(code=unique(res$code), decision = 0)
table(codes$decision)

# (1) Non-significant wald ratio
temp <- subset(res, nsnp == 1 & pval > threshold)$code
codes$decision[codes$code %in% temp] <- 1
table(codes$decision)


# (2) Significant wald ratio
temp <- subset(res, nsnp == 1 & pval <= threshold)$code
codes$decision[codes$code %in% temp] <- 2
table(codes$decision)

# (3) ivw not sig

temp <- subset(res, method == "Inverse variance weighted" & pval > threshold)$code
codes$decision[codes$code %in% temp] <- 3
table(codes$decision)

### Followup

# (4) 2-9 inst, no heterogeneity
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 2 & nsnp < 10)$code
het_followup <- subset(het, code %in% ivw_followup)
temp <- subset(het_followup, Q_pval > 0.05/nrow(het_followup))$code
codes$decision[codes$code %in% temp] <- 4
table(codes$decision)

# (5) heterogeneity

temp <- subset(het_followup, Q_pval <= 0.05/nrow(het_followup))$code
codes$decision[codes$code %in% temp] <- 5
table(codes$decision)


# (6) 10-49 inst, no het or high sign
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 10 & nsnp < 50)
het_followup <- subset(het, code %in% ivw_followup$code)
temp1 <- subset(het_followup, Q_pval > 0.05/nrow(het_followup))$code
temp2 <- subset(res, code %in% ivw_followup$code & method == "Sign concordance test")
temp2 <- subset(temp2, pval < 0.05/nrow(temp2))$code
temp3 <- subset(ivw_followup, code %in% c(temp1, temp2))$code

codes$decision[codes$code %in% temp3] <- 6
table(codes$decision)



# (7) high het/ low sign and median is not sig
temp4 <- subset(ivw_followup, ! code %in% c(temp1, temp2))
temp5 <- subset(res, code %in% temp4$code & method == "Simple median")
temp6 <- subset(temp5, pval > 0.05/nrow(temp5))$code
codes$decision[codes$code %in% temp6] <- 7
table(codes$decision)


# (8) high het/ low sign and median is sig
temp7 <- subset(temp5, pval <= 0.05/nrow(temp5))$code
codes$decision[codes$code %in% temp7] <- 8
table(codes$decision)


# (9) 50+ inst, low sign or high het
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 50)

temp <- subset(res, code %in% ivw_followup$code & method == "Sign concordance test")
temp2 <- subset(het, code %in% ivw_followup$code)

temp3 <- subset(temp, pval < 0.05/nrow(temp))$code
temp4 <- subset(temp2, Q_pval > 0.05/nrow(temp2))$code

temp5 <- subset(ivw_followup, !code %in% c(temp3, temp4))$code
codes$decision[codes$code %in% temp5] <- 9
table(codes$decision)


# (10) 50+ inst, high sign or low het
temp6 <- subset(ivw_followup, code %in% c(temp3, temp4))$code
codes$decision[codes$code %in% temp6] <- 10
table(codes$decision)
codes$sig <- codes$decision %in% c(2,4,6,8,10)
codes$msig <- codes$decision %in% c(4,6,8,10)

save(codes, file="../results/mrbase_sig_codes.rdata")



load("../results/mrbase_sig_codes.rdata")
load("../results/mrbase_tophits_full.rdata")
res$code <- paste(res$id.exposure, res$id.outcome)
res <- inner_join(res, codes) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"))
res$trait <- strsplit(res$exposure, split="\\|") %>% sapply(function(x) x[1]) %>% gsub(" $", "", .)

ressig <- subset(res, sig)
ressig1 <- subset(res, decision == 2)
ressig1 %>% filter(grepl("smok", exposure, ignore.case=TRUE))
ressigm <- subset(res, msig)
ressigm %>% group_by(trait) %>% summarise(n=n())
ressigm %>% group_by(outcome) %>% summarise(n=n())
subset(res, decision == 2) %>% group_by(trait) %>% summarise(n=n())


a <- subset(res, nsnp > 1, select=-c(id.exposure, id.outcome, method, chunk, code, sig, exposure))

a <- subset(res, nsnp > 1, select=c(trait, outcome, nsnp, b, se, pval, decision, msig)) %>%
	arrange(desc(msig), pval)
write.csv(a, file="../results/trait-cpg-sig-codes.csv")

