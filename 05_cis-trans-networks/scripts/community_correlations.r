library(dplyr)
library(igraph)

getcor <- function(creg, tcpg)
{
	cor(betas[which(rownames(betas) == creg)[1],], betas[which(rownames(betas) == tcpg)[1],], use="pair")
}

load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj")

load("../results/creg_tcpg.rdata")
load("../results/graph.rdata")

dat2 <- group_by(dat, creg) %>%
	arrange(tcpg_pval) %>%
	filter(!duplicated(tcpg_chr))

cregs <- unique(dat2$creg)
tcpgs <- unique(dat2$tcpg)

index <- rownames(norm.beta.random) %in% unique(c(cregs, tcpgs))

betas <- norm.beta.random[index,]

corres <- group_by(dat2, creg, tcpg) %>%
	mutate(r=getcor(creg,tcpg))

# # Test
# dat2_perm <- dat2
# dat2_perm$tcpg <- sample(dat2_perm$tcpg)
# corres_perm <- group_by(dat2_perm, creg, tcpg) %>%
# 	mutate(r=getcor(creg,tcpg))
# sum(is.na(corres_perm$r))
# mean(corres_perm$r^2, na.rm=TRUE)
# mean(corres$r^2, na.rm=TRUE)

# summary(lm(r ~ waldratio, corres))
# summary(lm(r ~ waldratio, corres_perm))

nperm <- 100
mean_cor_perms <- array(0, nperm)
for(i in 1:nperm)
{
	message(i)
	dat2_perm$tcpg <- sample(dat2_perm$tcpg)
	corres_perm <- group_by(dat2_perm, creg, tcpg) %>%
		mutate(r=getcor(creg,tcpg))
	mean_cor_perms[i] <- mean(corres_perm$r^2, na.rm=TRUE)
	message(mean_cor_perms[i])
}

save(corres, mean_cor_perms, file="../results/corres.rdata")



