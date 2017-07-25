arguments<-commandArgs(T)
no<-as.numeric(arguments[1])


load("../results/enrichments/snpcontrolsets.rdata")
library(ggplot2)

labs <- data.frame(table(f.all$groups))
m<-match(f.all$groups,labs[,1])
f.all<-data.frame(f.all,Ncatgroups=as.character(labs[m,-1]))

f.all$groups<-as.factor(f.all$groups)
f.all$groups <- factor(f.all$groups, levels = f.all$groups[order(as.numeric(as.character(f.all$Ncatgroups)),decreasing=T)])

f.all$MAF<-as.numeric(as.character(f.all$MAF))
f.all$nproxies<-as.numeric(as.character(f.all$nproxies))
f.all$tssdist<-as.numeric(as.character(f.all$tssdist))
f.all$closest450kdistance<-as.numeric(as.character(f.all$closest450kdistance))

quantile(f.all$MAF)
#     0%     25%     50%     75%    100% 
#0.01000 0.04175 0.13520 0.29620 0.50000 
quantile(f.all$nproxies)
#  0%  25%  50%  75% 100% 
#   0    1   15   47 1773 
quantile(f.all$tssdist)
#       0%       25%       50%       75%      100% 
#        0     11336     46753    211620 154984156 
quantile(f.all$closest450kdistance)
#     0%     25%     50%     75%    100% 
#      0    2582    8916   27030 2654985 
 
#applyquintiles <- function(x) {
# cut(x, breaks=c(quantile(x, probs = seq(0, 1, by = 0.20))), 
#     labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
#}
#f.all$MAFquantile <- sapply(f.all$MAF, applyquintiles(x=f.all$MAF))
#table(df$Quintile)


f.all$MAFquantile <- cut(f.all$MAF, breaks=c(quantile(f.all$MAF,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$MAFquantile)

f.all$nproxiesquantile <- cut(f.all$nproxies, breaks=c(quantile(f.all$nproxies,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)
table(f.all$nproxiesquantile)

f.all$tssdistquantile <- cut(f.all$tssdist, breaks=c(quantile(f.all$tssdist,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$tssdistquantile)

f.all$closest450kdistancequantile <- cut(f.all$closest450kdistance, breaks=c(quantile(f.all$closest450kdistance,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$closest450kdistancequantile)

f.all$groups_quantiles<-paste(f.all$MAFquantile,f.all$nproxiesquantile,f.all$tssdistquantile,f.all$closest450kdistancequantile)

##select groups with mQTLs
t<-table(f.all$groups_quantiles,f.all$mQTL)
m<-match(f.all$groups_quantiles,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols=t2[,1], Nmqtl=t2[,2])
t<-t[t[,2]>0,]

dim(t)
#625

#controls<-data.frame(SNP=f.all$SNP)

controllist<-list()
mqtllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
f.all.subgroup<-f.all[which(f.all$groups_quantiles==row.names(t)[i]&f.all$mQTL==FALSE),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$mQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist",no,".rdata",sep=""))

#w<-which(controls$SNP%in%id)
#w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$mQTL==TRUE)
#snp<-rep(NA,length(f.all$SNP))
#snp[w]<-0
#snp[w1]<-1
#controls<-data.frame(controls,SNP=snp)
#}

#save(controls,file=paste("../results/enrichments/controls",no,".rdata",sep=""))

#####

###
#cisSNPs


t<-table(f.all$groups_quantiles,f.all$cismQTL)
m<-match(f.all$groups_quantiles,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Nciscontrols=t2[,1], Ncismqtl=t2[,2])
t<-t[t[,2]>0,]

dim(t)

#controls<-data.frame(SNP=f.all$SNP)
controllist<-list()
mqtllist<-list()


for (i in 1:dim(t)[1]){
cat(i,"\n")
f.all.subgroup<-f.all[which(f.all$groups_quantiles==row.names(t)[i]&f.all$cismQTL==TRUE),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$mQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}
save(mqtllist,controllist,file=paste("../results/enrichments/controlslist_cis",no,".rdata",sep=""))


#w<-which(controls$SNP%in%id)
#w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$cismQTL==TRUE)

#snp<-rep(NA,length(f.all$SNP))
#snp[w]<-0
#snp[w1]<-1
#controls<-data.frame(controls,SNP=snp)
#}

#save(controls,file=paste("../results/enrichments/controls_cis",no,".rdata",sep=""))

#transSNPs

t<-table(f.all$groups_quantiles,f.all$transmQTL)
m<-match(f.all$groups_quantiles,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ntranscontrols=t2[,1], Ntransmqtl=t2[,2])
t<-t[t[,2]>0,]

dim(t)



#controls<-data.frame(SNP=f.all$SNP)
controllist<-list()
mqtllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
f.all.subgroup<-f.all[which(f.all$groups_quantiles==row.names(t)[i]&f.all$transmQTL==TRUE),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$mQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}
save(mqtllist,controllist,file=paste("../results/enrichments/controlslist_trans",no,".rdata",sep=""))


#w<-which(controls$SNP%in%id)
#w1<-which(f.all$groups_quantiles==row.names(t)[i]&f.all$transmQTL==TRUE)
#snp<-rep(NA,length(f.all$SNP))
#snp[w]<-0
#snp[w1]<-1
#controls<-data.frame(controls,SNP=snp)
#}

#save(controls,file=paste("../results/enrichments/controls_trans",no,".rdata",sep=""))

#save(f.all,file="../results/enrichments/snpcontrolsets.rdata")
