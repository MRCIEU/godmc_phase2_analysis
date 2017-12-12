arguments<-commandArgs(T)
no<-as.numeric(arguments[1])

library(ggplot2)
library(gridExtra)

load("../results/enrichments/snpcontrolsets_selection.rdata")

w<-which(is.na(f.all$snp_cis))
f.all$snp_cis<-as.character(f.all$snp_cis)

f.all$snp_cis[w]<-"no_mqtl"
table(f.all$snp_cis)

#ambivalent      FALSE    no_mqtl       TRUE 
#      2842      13516    9852470     216244 
f.all$snp_cis<-as.factor(f.all$snp_cis)

#remove snps
flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(f.all$SNP%in%flip[,1])
length(w) #2911
if(length(w)>0){
f.all<-f.all[-w,]}

indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
w<-which(f.all$SNP%in%indels[,1]) 
length(w) #2805
if(length(w)>0){
f.all<-f.all[-w,]}

f.all$MAF<-as.numeric(as.character(f.all$MAF))
f.all$nproxies<-as.numeric(as.character(f.all$nproxies))
f.all$tssdist<-as.numeric(as.character(f.all$tssdist))
f.all$closest450kdistance<-as.numeric(as.character(f.all$closest450kdistance))
f.all$GC_freq<-as.numeric(as.character(f.all$GC_freq))
f.all$CpG_freq<-as.numeric(as.character(f.all$CpG_freq))

table(f.all$mQTL)

#    FALSE    TRUE 
#9852470  232602 
 
dim(f.all[,1:3])
#[1] 10079356        3
dim(unique(f.all[,1:3]))
#3967661       3

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
quantile(f.all$GC_freq)
#       0%       25%       50%       75%      100% 
#0.0000000 0.3707614 0.3982461 0.4368927 1.0000000 
quantile(f.all$CpG_freq)
#         0%         25%         50%         75%        100% 
#0.000000000 0.006098446 0.008068602 0.011822609 1.000000000 

#applyquintiles <- function(x) {
# cut(x, breaks=c(quantile(x, probs = seq(0, 1, by = 0.20))), 
#     labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
#}
#f.all$MAFquantile <- sapply(f.all$MAF, applyquintiles(x=f.all$MAF))
#table(df$Quintile)



f.all$MAFquantile <- cut(log(f.all$MAF+0.001), breaks=c(quantile(f.all$MAF,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$MAFquantile)

f.all$nproxiesquantile <- cut(log(f.all$nproxies+0.001), breaks=c(quantile(f.all$nproxies,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)
table(f.all$nproxiesquantile)

f.all$tssdistquantile <- cut(log(abs(f.all$tssdist)+0.001), breaks=c(quantile(f.all$tssdist,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$tssdistquantile)

f.all$closest450kdistancequantile <- cut(log(abs(f.all$closest450kdistance)+0.001), breaks=c(quantile(f.all$closest450kdistance,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$closest450kdistancequantile)

f.all$CpG_freqquantile <- cut(f.all$CpG_freq, breaks=c(quantile(f.all$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$CpG_freqquantile)

f.all$GC_freqquantile <- cut(f.all$GC_freq, breaks=c(quantile(f.all$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$GC_freqquantile)

f.all$groups_quantiles<-paste(f.all$MAFquantile,f.all$nproxiesquantile,f.all$tssdistquantile,f.all$closest450kdistancequantile)
f.all$groups_quantiles2<-paste(f.all$MAFquantile,f.all$nproxiesquantile,f.all$tssdistquantile,f.all$closest450kdistancequantile,f.all$CpG_freqquantile,f.all$GC_freqquantile)

p1 <- ggplot(f.all, aes(MAF, fill=mQTL)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all, aes(nproxies, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all, aes(tssdist, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all,aes(closest450kdistance,fill=mQTL)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all, aes(CpG_freq, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all, aes(GC_freq, fill=mQTL)) + 
geom_density(alpha = 0.2)


pdf("../images/enrichmentpropertiesdensity.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

p1 <- ggplot(f.all, aes(log(MAF+0.001), fill=mQTL)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all, aes(log(nproxies+0.001), fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,10)
p3 <- ggplot(f.all, aes(log(abs(tssdist)+0.001), fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,15)
p4<-ggplot(f.all,aes(log(abs(closest450kdistance)+0.001),fill=mQTL)) +
geom_density(alpha = 0.2) +
xlim(0,15)
p5 <- ggplot(f.all, aes(CpG_freq, fill=mQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all, aes(GC_freq, fill=mQTL)) + 
geom_density(alpha = 0.2)


pdf("../images/enrichmentpropertiesdensitylogged.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

pdf("nproxies.pdf",height=6,width=6)
par(mfrow=c(2,2))
hist(log(f.all$nproxies+0.001),main="")
hist(log(f.all$nproxies+0.01),main="")
hist(log(f.all$nproxies+0.1),main="")
dev.off()

f.all.cis<-f.all[which(f.all$snp_cis!="FALSE"&f.all$snp_cis!="ambivalent"),]

p1 <- ggplot(f.all.cis, aes(MAF, fill=snp_cis)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all.cis, aes(nproxies, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all.cis, aes(tssdist, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all.cis,aes(closest450kdistance,fill=snp_cis)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all.cis, aes(CpG_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all.cis, aes(GC_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensitycis.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

f.all.trans<-f.all[which(f.all$snp_cis!="TRUE"&f.all$snp_cis!="ambivalent"),]
p1 <- ggplot(f.all.trans, aes(MAF, fill=snp_cis)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all.trans, aes(nproxies, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all.trans, aes(tssdist, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all.trans,aes(closest450kdistance,fill=snp_cis)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all.trans, aes(CpG_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all.trans, aes(GC_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensitytrans.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

f.all.amb<-f.all[which(f.all$snp_cis!="FALSE"&f.all$snp_cis!="TRUE"),]
p1 <- ggplot(f.all.amb, aes(MAF, fill=snp_cis)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all.amb, aes(nproxies, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all.amb, aes(tssdist, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all.amb,aes(closest450kdistance,fill=snp_cis)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all.amb, aes(CpG_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all.amb, aes(GC_freq, fill=snp_cis)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensityambivalent.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

##select groups with mQTLs

t<-table(f.all$groups_quantiles2,f.all$mQTL)
dim(t)
#11884
m<-match(f.all$groups_quantiles2,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_6prop=t2[,1], Nmqtl_6prop=t2[,2])
t<-data.frame(t[t[,2]>0,])
dim(t)
#[1] 9088    2

t$req10<-10*as.numeric(t[,2])
w<-which(t[,1]-t[,3]<0)
cat(nrow(t[w,]),"controlgroups have less than 10 times mqtl snps","\n")

t$req5<-5*as.numeric(t[,2])
w<-which(t[,1]-t[,4]<0)
cat(nrow(t[w,]),"controlgroups have less than 5 times mqtl snps","\n")

w<-which(t[,1]-t[,2]<0)
cat(nrow(t[w,]),"controlgroups < mqtl snps","\n")

df1<-data.frame(no=t[,1],type="nonmqtls")
df2<-data.frame(no=t[,2],type="mqtls")
df<-rbind(df1,df2)
df$group<-cut(df$no, breaks=c(seq(0, 5000, by = 50)))
df<-df[which(df$no<20000),]
df2<-data.frame(table(df$type,df$group))

df3<-df2[which(df2$Var1=="mqtls"),]
p1 <- ggplot(df3, aes(x=as.factor(df3$Var2),y=Freq)) +
geom_bar(stat="identity") +
theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90, vjust=0.5, size=4)) +
xlab("groupsizes")
ggsave(plot=p1, file="groupsizesmqtls_6properties.pdf", width=7, height=7)

df3<-df2[which(df2$Var1=="nonmqtls"),]
p1 <- ggplot(df3, aes(x=as.factor(df3$Var2),y=Freq)) +
geom_bar(stat="identity") +
theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90, vjust=0.5, size=4)) +
xlab("groupsizes")
ggsave(plot=p1, file="groupsizesnonmqtls_6properties.pdf", width=7, height=7)



controllist<-list()
mqtllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[,1]>5&t[,2]>0){
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$mQTL==FALSE&f.all$snpchr!="chrX"),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$mQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop",no,".rdata",sep=""))



#####

###
#cisSNPs


t<-table(f.all$groups_quantiles2,f.all$snp_cis=="TRUE")
dim(t)
m<-match(f.all$groups_quantiles2,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_6prop_cis=t2[,1], Nmqtl_6prop_cis=t2[,2])
t<-data.frame(t[t[,2]>0,])
dim(t)


controllist<-list()
mqtllist<-list()


for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[,1]>5&t[,2]>0){
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis!="TRUE"),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis=="TRUE")
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop_cis",no,".rdata",sep=""))


#transSNPs

t<-table(f.all$groups_quantiles2,f.all$snp_cis=="FALSE")
dim(t)
m<-match(f.all$groups_quantiles2,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_6prop_trans=t2[,1], Nmqtl_6prop_trans=t2[,2])
t<-data.frame(t[t[,2]>0,])
dim(t)


controllist<-list()
mqtllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[,1]>5&t[,2]>0){
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis!="FALSE"),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis=="FALSE")
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop_trans",no,".rdata",sep=""))

#ambivalentSNPs

t<-table(f.all$groups_quantiles2,f.all$snp_cis=="ambivalent")
dim(t)
m<-match(f.all$groups_quantiles2,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_6prop_amb=t2[,1], Nmqtl_6prop_amb=t2[,2])
t<-data.frame(t[t[,2]>0,])
dim(t)


controllist<-list()
mqtllist<-list()

for (i in 1:dim(t)[1]){
cat(i,"\n")
if(t[,1]>5&t[,2]>0){
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis!="ambivalent"),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$snp_cis=="ambivalent")
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop_ambivalent",no,".rdata",sep=""))


