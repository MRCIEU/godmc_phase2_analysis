library(dplyr)
load("results/16/16_conditional.rdata")

group_by(conditional, cpg, cis) %>% summarise(n=n()) %>% arrange(desc(n))

table(conditional$cis)
sum(!conditional$cis) / nrow(conditional)

temp <- filter(conditional, cis) %>% group_by(cpg, cpgchr=Chr) %>% summarise()
conditional <- left_join(conditional, temp)
conditional$cischr <- conditional$Chr == conditional$cpgchr
conditional$cischr[is.na(conditional$cischr)] <- FALSE

temp2 <- group_by(conditional, cpg, cischr) %>% summarise(n=n())

sum(!temp2$cischr) / nrow(temp2)


a <- bind_rows(
	filter(conditional, cis, pJ < 1e-8),
	filter(conditional, !cis, pJ < 1e-14)
)
dim(a)

b <- a %>% group_by(cpg) %>% summarise(n=n())

mean(b$n)
median(b$n)
table(b$n) %>% sort(descending=TRUE) %>% sum
summary(b$n)



b <- filter(a, cis) %>% group_by(cpg) %>% summarise(n=n(), d=max(bp) - min(bp))

load("results/16/16_clumped.rdata")
clmin <- filter(clumped, cis) %>% group_by(cpg) %>% summarise(minp=min(pval))

b <- inner_join(b, clmin)
b$minp <- pmax(1e-300, b$minp)

cor(b$n, b$d)
cor(b$n, -log10(b$minp))


