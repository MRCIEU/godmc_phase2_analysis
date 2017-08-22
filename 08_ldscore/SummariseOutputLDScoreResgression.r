dat<-read.table("LDScore_MergedOutput.txt", header = TRUE, row.names = NULL, sep = "\t")

## calc CI around h2 estimate

upper<-dat$h2+1.96*dat$se
lower<-dat$h2-1.96*dat$se

dat<-cbind(dat[,1:6],lower,upper,dat[,-c(1:6)])

## indicator which probes have nonzero estimates of h2
nonZero<-rep(0, nrow(dat))
nonZero[which(lower > 0)]<-1

## do probes with nozero h2 have higher estimates in the van Dongen paper?
pdf("BoxplotVanDongenEstimatesSplitBynonZeroLDSR.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
boxplot(dat$h2_total ~ nonZero, ylab = "h2 total", xlab = "nonZero LD score estimate", names = c("No", "Yes"))
mtext(paste("P = ", signif(wilcox.test(dat$h2_total ~ nonZero)$p.value,3), sep = ""), side = 3, adj = 1)
boxplot(dat$h2_SNPs ~ nonZero, ylab = "h2 SNPs", xlab = "nonZero LD score estimate", names = c("No", "Yes"))
mtext(paste("P = ", signif(wilcox.test(dat$h2_SNPs ~ nonZero)$p.value,3), sep = ""), side = 3, adj = 1)
dev.off()

## plot against estimates from van Dongen
pdf("ScatterplotplotVanDongenEstimatesLDSR.pdf", width = 10, height = 5)

par(mfrow = c(1,2))

plot(dat$h2, dat$h2_total, xlab = "LD score estimate", ylab = "h2 total", pch = 16, col = factor(nonZero))
plot(dat$h2, dat$h2_SNPs, xlab = "LD score estimate", ylab = "h2 SNPs", pch = 16, col = factor(nonZero))
dev.off()
