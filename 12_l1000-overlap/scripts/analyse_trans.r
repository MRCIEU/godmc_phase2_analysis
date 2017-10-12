l <- list()
for(i in 1:10)
{
	load(paste0("../results/perturbations_Zscores_perm", i, ".Rdata"))
	l[[i]] <- output
}

load("../results/perturbations_Zscores.Rdata")

mean(abs(output$Z))
sapply(l, function(x) mean(abs(x$Z)))

sum(abs(output$Z) > qnorm(1e-5, low=FALSE))
sapply(l, function(x) sum(abs(x$Z) > qnorm(1e-5, low=FALSE)))

max(abs(output$Z))
sapply(l, function(x) max(abs(x$Z)))


library(tidyr)

org <- function(fn)
{
	load(fn)
	require(tidyr)
	require(dplyr)
	output <- separate(output, sig_info, c("code", "pgene", "id"), sep=":", rem=FALSE) %>%
		group_by(cpg, snp, snp_hgnc_symbol, cpg_hgnc_symbol, code, pgene) %>%
		summarise(Z = mean(Z), n=n())
}


l <- list()
for(i in 1:10)
{
	message(i)
	output <- org(paste0("../results/perturbations_Zscores_perm", i, ".Rdata"))
	l[[i]] <- output
	l[[i]]$perm <- i
}

l <- bind_rows(l)

output <- org("../results/perturbations_Zscores.Rdata")

x <- subset(output, code == code[1])

ttest <- group_by(output, code) %>%
	do({
		x <- .
		message(x$code[1])
		perm <- subset(l, code == x$code[1])$Z %>% abs
		res <- t.test(abs(x$Z), perm)
		return(data_frame(statistic=res$statistic, df=res$parameter, pval=res$p.value))
	})

m <- group_by(l, code)



library(dplyr)
library(tidyr)
load("../results/same_signature/perturbations_Zscores_pairwise.Rdata")


org <- function(fn)
{
	load(fn)
	require(tidyr)
	require(dplyr)
	output <- separate(output, sig_info, c("code", "id"), sep=":", rem=FALSE)
}


l <- list()
for(i in 1:10)
{
	message(i)
	output <- org(paste0("../results/same_signature/perturbations_Zscores_pairwise_perm", i, ".Rdata"))
	l[[i]] <- output
	l[[i]]$perm <- i
}

output <- org("../results/same_signature/perturbations_Zscores_pairwise.Rdata")


ap <- function(fn)
{
	data_frame(val=c(fn(output), sapply(l, fn)), what=c("TRUE", rep("Perm", length(l)))) %>%
	arrange(desc(val))
}

ap(nrow)
ap(function(x) mean(abs(x$Z_SNP)))
ap(function(x) mean(abs(x$Z_CpG)))


temp <- group_by(output, id, code) %>%
summarise(n=n(), what="True")

temp2 <- lapply(l, function(x)
{
	group_by(x, id, code) %>%
	summarise(n=n(), what=first(perm))
}) %>% bind_rows

table(temp2$what)
nrow(temp)

group_by(temp2, what) %>%
summarise(min=min(n), max=max(n), median=median(n), mean(n), sd(n))

summary(temp$n)
sd(temp$n)
