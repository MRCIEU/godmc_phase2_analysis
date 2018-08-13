library(TwoSampleMR)
library(tidyverse)


load("../data/snps_gwas.rdata")
## loads a
gkeep <- filter(a, grepl("mr", data_source.exposure))
gkeeps <- gkeep %>% group_by(id.exposure, trait, exposure) %>% summarise(n=n()) %>% arrange(desc(n)) %>% ungroup %>% filter(!duplicated(trait))
rm(a)


# Filtering 06 results


load("../results/mrbase_sig.rdata")
## loads sig
sig$code <- paste(sig$id.exposure, sig$id.outcome)

# 1. wald ratios
code_wr <- subset(sig, nsnp == 1 & pval < 1.4e-7 & id.exposure %in% gkeeps$id.exposure)$code



# 2. 2-25 SNPs
## if low Q stat, use IVW

load("../results/mrbase_sig_mhc.rdata")
## loads res

mrbase_sig_mhc <- subset(res, id.exposure %in% gkeeps$id.exposure & pval < 1.4e-7)
rm(res)
mrbase_sig_mhc <- subset(mrbase_sig_mhc, what2 == "all" & nsnp > 1 & nsnp < 25 & is.na(mrbase_sig_mhc$what))
mrbase_sig_mhc$code <- paste(mrbase_sig_mhc$id.exposure, mrbase_sig_mhc$id.outcome)

code_lowq_25 <- subset(mrbase_sig_mhc, Q_pval > 0.05/nrow(mrbase_sig_mhc))$code


## if high Q stat, use mode
load("../results/tophits_followup.rdata")
## loads res
tophits_followup <- res

temp <- select(gkeeps, exposure, id.exposure)
query_codes <- subset(mrbase_sig_mhc, Q_pval < 0.05/nrow(mrbase_sig_mhc))$code

out <- inner_join(tophits_followup, temp) %>% mutate(id.outcome = outcome, code = paste(id.exposure, id.outcome)) %>% filter(code %in% query_codes & method == "Simple mode" & pval < 1.4e-7)
code_hiq_25 <- out$code
rm(res, temp, query_codes)

# 3. 50+ SNPs
## check if sct is significant

load("../../06_mr-gwas-cpg/results/mrbase_sig_mhc_sign.rdata")

mult_code <- subset(res, nsnp >= 50 & what2 == "all" & method == "Inverse variance weighted" & pval < 1.4e-7)$code
code_50 <- subset(res, code %in% mult_code & what2 == "all" & method == "Sign concordance test" & pval < 0.05/length(mult_code))$code

mult <- group_by(mult, code, exposure) %>%
	summarise(same_sign=sign(b[1]) == sign(b[2]), sig1=pval[2] < 0.05/nrow(mult)*2, sig2 = pval[2] < 0.05, sig3 = pval[1] < threshold2)
table(mult$same_sign, mult$sig2, mult$sig3)
group_by(mult, exposure) %>% summarise(n=n(), nsig1=sum(sig1), nsig2=sum(sig2), nsig3=sum(sig2 & sig3)) %>% arrange( nsig2) %>% as.data.frame
sum(mult$sig2 & mult$sig3)

ressig1 <- subset(res, (method == "Wald ratio" | method == "Inverse variance weighted") & what2=="all" & pval < threshold2)
ressig1$sig <- ressig1$code %in% subset(mult, sig2 & sig3)$code
ressig1 <- inner_join(ressig1, cpgpos, by=c("outcome"="cpg"))
ressig1$chr <- as.numeric(gsub("chr", "", ressig1$cpgchr))
ressig1$trait <- strsplit(ressig1$exposure, split="\\|") %>% sapply(function(x) x[1]) %>% gsub(" $", "", .)
# ressig2 <- subset(res, method == "Inverse variance weighted" & pval < threshold2)
# ressig2$nom <- strsplit(ressig2$exposure, split=" ") %>% sapply(function(x) x[1])
ressig1$pval[ressig1$pval < 1e-100] <- 1e-100



