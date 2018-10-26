library(data.table)
load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped<-clumped[which(clumped$cis==TRUE & clumped$pval<1e-8 | clumped$cis==FALSE & clumped$pval<1e-14),]
table(clumped$cis)
# FALSE   TRUE 
# 23117 248607 
length(unique(clumped$cpg))
#190102

trans_interchrom<-clumped[which(clumped$snpchr!=clumped$cpgchr),]
#18584
nrow(trans_interchrom)

res <- rep(0,1000)
for(i in 1:1000)
{
	message(i)
	load(paste0("~/repo/godmc_phase2_analysis/data/hi-c/permutations/nodups_data_perm_",i, ".Rdata"))
	res[i] <- length(unique(nodups_data$code))
}
perm<-median(res)
#473
max(res) #547
load(paste0("~/repo/godmc_phase2_analysis/data/hi-c/nodups.data.Rdata"))
real<-length(unique(nodups_data$code))
#637
sort(res)

#hic<-matrix(c(1175,517,18644-1175,18644-517),nrow=2,dimnames=list(c("real","permuted"),c("overlap","no overlap")))

hic<-matrix(c(real,perm,nrow(trans_interchrom)-real,nrow(trans_interchrom)-perm),nrow=2,dimnames=list(c("real","permuted"),c("overlap","no overlap")))

fisher.test(hic)


#	Fisher's Exact Test for Count Data

#data:  hic
#p-value = 6.506e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.202451 1.536728
#sample estimates:
#odds ratio 
#  1.359012 

data=as.data.table(clumped)
data[,cpg_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("cpgchr","cpgpos")]
hic_cpg<-unique(nodups_data$CpG)
cl<-data[which(data$cpg%in%hic_cpg),]
df<-unique(data.frame(cl$cpg,cl$cpg_cis))
table(df[,2])

#ambivalent      FALSE 
#       368        234





