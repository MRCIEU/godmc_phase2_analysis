res <- rep(0,1000)
for(i in 1:1000)
{
	message(i)
	load(paste0("permutations/nodups_data_perm_",i, ".Rdata"))
	res[i] <- length(unique(nodups_data$code))
}
perm<-median(res)
#473
load(paste0("nodups.data.Rdata"))
real<-length(unique(nodups_data$code))
#637
sort(res)

#hic<-matrix(c(1175,517,18644-1175,18644-517),nrow=2,dimnames=list(c("real","permuted"),c("overlap","no overlap")))

hic<-matrix(c(real,perm,18644-1175,18644-517),nrow=2,dimnames=list(c("real","permuted"),c("overlap","no overlap")))

fisher.test(hic)

#	Fisher's Exact Test for Count Data

#data:  hic
#p-value = 5.426e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.236446 1.580211
#sample estimates:
#odds ratio 
#  1.397429 