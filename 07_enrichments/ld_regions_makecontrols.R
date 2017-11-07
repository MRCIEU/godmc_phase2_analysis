arguments<-commandArgs(T)
no<-as.numeric(arguments[1])

library(ggplot2)
library(gridExtra)

load("../results/enrichments/snpcontrolsetsGC_CpGcontent.rdata")
f.all<-r.all


#labs <- data.frame(table(f.all$groups))
#m<-match(f.all$groups,labs[,1])
#f.all<-data.frame(f.all,Ncatgroups=as.character(labs[m,-1]))

#f.all$groups<-as.factor(f.all$groups)
#f.all$groups <- factor(f.all$groups, levels = f.all$groups[order(as.numeric(as.character(f.all$Ncatgroups)),decreasing=T)])
flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(f.all$SNP%in%flip[,1])
f.all<-f.all[-w,]

f.all$MAF<-as.numeric(as.character(f.all$MAF))
f.all$nproxies<-as.numeric(as.character(f.all$nproxies))
f.all$tssdist<-as.numeric(as.character(f.all$tssdist))
f.all$closest450kdistance<-as.numeric(as.character(f.all$closest450kdistance))
f.all$GC_freq<-as.numeric(as.character(f.all$GC_freq))
f.all$CpG_freq<-as.numeric(as.character(f.all$CpG_freq))

#f.all$region<-paste(f.all$snpchr,f.all$min,f.all$max,sep="_")
#o<-order(f.all$mQTL,decreasing=T)
#f.all<-f.all[o,]
#region<-unique(f.all$region)
#m<-match(region,f.all$region)
#f.all2<-f.all[m,]

table(f.all$mQTL)

#  FALSE    TRUE 
#9849491  232670 

#table(f.all2$mQTL)

#  FALSE    TRUE 
#3788302  180168 

dim(f.all[,1:3])
#10082161
dim(unique(f.all[,1:3]))
#3968470

#dim(f.all2[,1:3])
#3968470
#dim(unique(f.all2[,1:3]))
#3968470

#f.all<-f.all2


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



f.all$MAFquantile <- cut(f.all$MAF, breaks=c(quantile(f.all$MAF,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$MAFquantile)

f.all$nproxiesquantile <- cut(f.all$nproxies, breaks=c(quantile(f.all$nproxies,probs = seq(0, 1, by = 0.25))), labels=c("0-25","25-50","50-75","75-100"), include.lowest=TRUE)
table(f.all$nproxiesquantile)

f.all$tssdistquantile <- cut(f.all$tssdist, breaks=c(quantile(f.all$tssdist,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$tssdistquantile)

f.all$closest450kdistancequantile <- cut(f.all$closest450kdistance, breaks=c(quantile(f.all$closest450kdistance,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$closest450kdistancequantile)

f.all$CpG_freqquantile <- cut(f.all$CpG_freq, breaks=c(quantile(f.all$CpG_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$CpG_freq)

f.all$GC_freqquantile <- cut(f.all$GC_freq, breaks=c(quantile(f.all$GC_freq,probs = seq(0, 1, by = 0.20))), labels=c("0-20","20-40","40-60","60-80","80-100"), include.lowest=TRUE)
table(f.all$GC_freq)

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

p1 <- ggplot(f.all, aes(MAF, fill=cismQTL)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all, aes(nproxies, fill=cismQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all, aes(tssdist, fill=cismQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all,aes(closest450kdistance,fill=cismQTL)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all, aes(CpG_freq, fill=cismQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all, aes(GC_freq, fill=cismQTL)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensitycis.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

p1 <- ggplot(f.all, aes(MAF, fill=transmQTL)) + 
geom_density(alpha = 0.2) 
p2 <- ggplot(f.all, aes(nproxies, fill=transmQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,250)
p3 <- ggplot(f.all, aes(tssdist, fill=transmQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,400000)
p4<-ggplot(f.all,aes(closest450kdistance,fill=transmQTL)) +
geom_density(alpha = 0.2) +
xlim(0,100000)
p5 <- ggplot(f.all, aes(CpG_freq, fill=transmQTL)) + 
geom_density(alpha = 0.2) +
xlim(0,0.25)
p6 <- ggplot(f.all, aes(GC_freq, fill=transmQTL)) + 
geom_density(alpha = 0.2)

pdf("../images/enrichmentpropertiesdensitytrans.pdf", width=7, height=7)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
dev.off()

##select groups with mQTLs
t<-table(f.all$groups_quantiles,f.all$mQTL)
m<-match(f.all$groups_quantiles,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_4prop=t2[,1], Nmqtl_4prop=t2[,2])
t<-data.frame(t[t[,2]>0,])

t$req10<-10*as.numeric(t[,2])
w<-which(t[,1]-t[,3]<0)
t[w,]

dim(t)
#[1] 9184    2

t$req5<-5*as.numeric(t[,2])
w<-which(t[,1]-t[,4]<0)
t[w,]

w<-which(t[,1]-t[,2]<0)
t[w,]

df1<-data.frame(no=t[,1],type="nonmqtls")
df2<-data.frame(no=t[,2],type="mqtls")
df<-rbind(df1,df2)
df$group<-cut(df$no, breaks=c(seq(0, 100000, by = 1000)))
df2<-data.frame(table(df$type,df$group))

df3<-df2[which(df2$Var1=="mqtls"),]
p1 <- ggplot(df3, aes(x=as.factor(df3$Var2),y=Freq)) +
geom_bar(stat="identity") +
theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90, vjust=0.5, size=6)) +
xlab("groupsizes")
ggsave(plot=p1, file="groupsizesmqtls_4properties.pdf", width=7, height=7)

df3<-df2[which(df2$Var1=="nonmqtls"),]
p1 <- ggplot(df3, aes(x=as.factor(df3$Var2),y=Freq)) +
geom_bar(stat="identity") +
theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90, vjust=0.5, size=6)) +
xlab("groupsizes")
ggsave(plot=p1, file="groupsizesnonmqtls_4properties.pdf", width=7, height=7)


dim(t)
#500

head(t)
                       
#                        FALSE TRUE
#  0-20 0-25 0-20 0-20   62640 1322
#  0-20 0-25 0-20 20-40  47610  638
#  0-20 0-25 0-20 40-60  37944  416
#  0-20 0-25 0-20 60-80  28438  203
#  0-20 0-25 0-20 80-100 15438   51
#  0-20 0-25 20-40 0-20  49858  951

t<-table(f.all$groups_quantiles2,f.all$mQTL)
dim(t)
m<-match(f.all$groups_quantiles2,row.names(t))
t2<-t[m,]
f.all<-data.frame(f.all,Ncontrols_6prop=t2[,1], Nmqtl_6prop=t2[,2])
t<-data.frame(t[t[,2]>0,])
dim(t)
#[1] 9184    2

t$req10<-10*as.numeric(t[,2])
w<-which(t[,1]-t[,3]<0)
t[w,]


t$req5<-5*as.numeric(t[,2])
w<-which(t[,1]-t[,4]<0)
t[w,]

w<-which(t[,1]-t[,2]<0)
t[w,]


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



#controllist<-list()
#mqtllist<-list()

#for (i in 1:dim(t)[1]){
#cat(i,"\n")
#if(t[,1]>5&t[,2]>0){
#f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$mQTL==FALSE),]
#id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
#controllist[[i]]<-id
#w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$mQTL==TRUE)
#mqtllist[[i]]<-f.all$SNP[w1]
#}}

#save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop",no,".rdata",sep=""))



#####

###
#cisSNPs


t<-table(f.all$groups_quantiles2,f.all$cismQTL)
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
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$cismQTL==FALSE),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$cismQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop_cis",no,".rdata",sep=""))


#transSNPs

t<-table(f.all$groups_quantiles2,f.all$transmQTL)
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
f.all.subgroup<-f.all[which(f.all$groups_quantiles2==row.names(t)[i]&f.all$transmQTL==FALSE),]
id<-f.all.subgroup[sample(nrow(f.all.subgroup), size=t[i,2], replace=FALSE),"SNP"]
controllist[[i]]<-id
w1<-which(f.all$groups_quantiles2==row.names(t)[i]&f.all$transmQTL==TRUE)
mqtllist[[i]]<-f.all$SNP[w1]
}}

save(mqtllist,controllist,file=paste("../results/enrichments/controlslist6prop_trans",no,".rdata",sep=""))



