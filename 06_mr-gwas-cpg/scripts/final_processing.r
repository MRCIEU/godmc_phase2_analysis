library(TwoSampleMR)
library(tidyverse)

load("../results/mrbase_tophits_full.rdata")

threshold <- 1.4e-7
res$code <- paste(res$id.exposure, res$id.outcome)
het$code <- paste(het$id.exposure, het$id.outcome)
het <- subset(het, method == "Inverse variance weighted")
plei$code <- paste(plei$id.exposure, plei$id.outcome)

codes <- data_frame(code=unique(res$code), outcome = 0)
table(codes$outcome)

# (1) Non-significant wald ratio
temp <- subset(res, nsnp == 1 & pval > threshold)$code
codes$outcome[codes$code %in% temp] <- 1
table(codes$outcome)


# (2) Significant wald ratio
temp <- subset(res, nsnp == 1 & pval <= threshold)$code
codes$outcome[codes$code %in% temp] <- 2
table(codes$outcome)

# (3) ivw not sig

temp <- subset(res, method == "Inverse variance weighted" & pval > threshold)$code
codes$outcome[codes$code %in% temp] <- 3
table(codes$outcome)

### Followup

# (4) 2-9 inst, no heterogeneity
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 2 & nsnp < 10)$code
het_followup <- subset(het, code %in% ivw_followup)
temp <- subset(het_followup, Q_pval > 0.05/nrow(het_followup))$code
codes$outcome[codes$code %in% temp] <- 4
table(codes$outcome)

# (5) heterogeneity

temp <- subset(het_followup, Q_pval <= 0.05/nrow(het_followup))$code
codes$outcome[codes$code %in% temp] <- 5
table(codes$outcome)


# (6) 10-49 inst, no het
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 10 & nsnp < 50)
het_followup <- subset(het, code %in% ivw_followup$code)
temp1 <- subset(het_followup, Q_pval > 0.05/nrow(het_followup))$code
temp2 <- subset(res, code %in% ivw_followup$code & method == "Sign concordance test")
temp2 <- subset(temp2, pval < 0.05/nrow(temp2))$code
temp3 <- subset(ivw_followup, code %in% c(temp1, temp2))$code

codes$outcome[codes$code %in% temp3] <- 6
table(codes$outcome)



# (7) high het/ low sign and median is not sig
temp4 <- subset(ivw_followup, ! code %in% c(temp1, temp2))
temp5 <- subset(res, code %in% temp4$code & method == "Simple median")
temp6 <- subset(temp5, pval > 0.05/nrow(temp5))$code
codes$outcome[codes$code %in% temp6] <- 7
table(codes$outcome)


# (8) high het/ low sign and median is sig
temp7 <- subset(temp5, pval <= 0.05/nrow(temp5))$code
codes$outcome[codes$code %in% temp7] <- 8
table(codes$outcome)


# (9) 50+ inst, low sign or high het
ivw_followup <- subset(res, method == "Inverse variance weighted" & pval <= threshold & nsnp >= 50)

temp <- subset(res, code %in% ivw_followup$code & method == "Sign concordance test")
temp2 <- subset(het, code %in% ivw_followup$code)

temp3 <- subset(temp, pval < 0.05/nrow(temp))$code
temp4 <- subset(temp2, Q_pval > 0.05/nrow(temp2))$code

temp5 <- subset(ivw_followup, !code %in% c(temp3, temp4))$code
codes$outcome[codes$code %in% temp5] <- 9
table(codes$outcome)


# (10) 50+ inst, high sign or low het
temp6 <- subset(ivw_followup, code %in% c(temp3, temp4))$code
codes$outcome[codes$code %in% temp6] <- 10
table(codes$outcome)
codes$sig <- codes$outcome %in% c(2,4,6,8,10)
codes$msig <- codes$outcome %in% c(4,6,8,10)

save(codes, file="../results/mrbase_sig_codes.rdata")
