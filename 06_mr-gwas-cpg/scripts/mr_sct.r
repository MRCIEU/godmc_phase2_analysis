library(TwoSampleMR)
library(dplyr)
library(ggplot2)

a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)
dat <- harmonise_data(a, b) %>% filter(mr_keep)

mr(dat)
summary(lm(beta.outcome ~ beta.exposure, dat))


set.seed(1)
dat2 <- dat
dat2$beta.outcome <- sample(dat2$beta.outcome)
out <- list()
for(i in 1:nrow(dat))
{
	message(i)
	dat3 <- dat
	dat3$se.outcome[i] <- dat3$se.outcome[i] / 5
	dat4 <- dat2
	dat4$se.outcome[i] <- dat4$se.outcome[i] / 5
	dat5 <- dat
	dat5$se.outcome[i] <- dat5$se.outcome[i] / 10
	dat6 <- dat2
	dat6$se.outcome[i] <- dat6$se.outcome[i] / 10

	res3 <- mr(dat3, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))
	res4 <- mr(dat4, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))
	res5 <- mr(dat5, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))
	res6 <- mr(dat6, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))
	res3$what <- res5$what <- "outlier"
	res4$what <- res6$what <- "null"
	res3$eff <- res4$eff <- 5
	res5$eff <- res6$eff <- 10
	out[[i]] <- rbind(res3,res4,res5,res6)
	out[[i]]$snp <- i
}
res <- bind_rows(out)




temp <- mr(dat, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))
temp$what <- "real"
temp$eff <- 0
temp$snp <- NA
res <- rbind(temp, res)

res <- subset(res, method != "Unweighted regression")
ress <- group_by(res, method, what, eff) %>%
summarise(nsig=sum(pval < 0.05))


set.seed(1)

# Extract SNP effects from MR-Base
a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)
dat2 <- harmonise_data(a,b)

# Permute SNP-outcome effects
dat2$beta.outcome <- sample(dat2$beta.outcome)

# Get MR estimates
temp2 <- mr(dat2, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))

# Plot
o2 <- mr_forest_plot(mr_singlesnp(dat2, all_method = c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode")))[[1]] +
labs(x="MR effect size estimates") + theme(axis.text.y=element_text(size=5))

dat3 <- dat2
dat3$se.outcome[25] <- dat3$se.outcome[25] / 5
temp3 <- mr(dat3, metho=c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode", "mr_uwr"))

dat3$outlier <- FALSE
dat3$outlier[25] <- TRUE
dat3$waldratio <- dat3$beta.outcome / dat3$beta.exposure

o3 <- mr_forest_plot(mr_singlesnp(dat3, all_method = c("mr_ivw_radial", "mr_sign", "mr_weighted_median", "mr_weighted_mode")))[[1]] + geom_point(data=dat3[dat3$outlier,], aes(y=SNP, x=waldratio), colour="black", size=2) + geom_point(data=dat3[dat3$outlier,], aes(y=SNP, x=waldratio), colour="yellow", size=1) +
labs(x="MR effect size estimates") + theme(axis.text.y=element_text(size=5))








save(res, dat, o3, o2, file="../results/mr_sct_bmi.rdata")


