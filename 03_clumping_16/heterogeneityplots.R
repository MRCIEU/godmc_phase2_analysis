#results are in /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16
library(tidyverse)

#962 chunks
#path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis"
#l<-list.files(path=paste0(path,"/results/16"),pattern=".txt.gz")

# Get cpg positions
#load(paste0(path,"/03_clumping_16/cpg_pos.rdata"))

#a14.cis.out<-data.frame()
#a14.trans.out<-data.frame()

#a8.cis.out<-data.frame()
#a8.trans.out<-data.frame()

#l=2
#for (i in 1:length(l)){

#cat(i,"\n")
#a <- read_tsv(paste0(path,"/results/16/16_", i, ".txt.gz"))
#a <- a %>% separate(MarkerName, into=c("snp", "cpg"), sep="_")
#a$snp2 <- a$snp
#a <- a %>% separate(snp2, into=c("snpchr", "snppos", "snptype"), sep=":")
#a$snppos <- as.numeric(a$snppos)
#a <- inner_join(a, cpgpos, by=c("cpg"))
#a$cis <- FALSE
#cis_radius <- 1000000
#a$cis[a$snpchr == a$cpgchr & (abs(a$snppos - a$cpgpos) <= cis_radius)] <- TRUE
#a<-data.frame(a)

#a14cis<-data.frame(a[which(a$cis==T & a$P.value<1e-14),"HetISq"],subset = "cis p<1e-14")
#a14trans<-data.frame(a[which(a$cis==F & a$P.value<1e-14),"HetISq"],subset = "trans p<1e-14")
#a8cis<-data.frame(a[which(a$cis==T & a$P.value<5e-8),"HetISq"],subset = "cis p<5e-8")
#a8trans<-data.frame(a[which(a$cis==F & a$P.value<5e-8),"HetISq"],subset = "trans p<5e-8")

#f<-function(x){
#y<-strsplit(x,split="")
#length(which(y[[1]]=="?"))}

#dir<-as.list(a$Direction)
#l3<-sapply(dir,f)

#a<-data.frame(a,l3)

#a14cis<-data.frame(a[which(a$cis==T & a$P.value<1e-14),],subset = "cis p<1e-14")
#a14trans<-data.frame(a[which(a$cis==F & a$P.value<1e-14),],subset = "trans p<1e-14")
#a8cis<-data.frame(a[which(a$cis==T & a$P.value<5e-8),],subset = "cis p<5e-8")
#a8trans<-data.frame(a[which(a$cis==F & a$P.value<5e-8),],subset = "trans p<5e-8")

#a14.cis.out<-rbind(a14.cis.out,a14cis)
#a8.cis.out<-rbind(a8.cis.out,a8cis)
#a14.trans.out<-rbind(a14.trans.out,a14trans)
#a8.trans.out<-rbind(a8.trans.out,a8trans)

#save(a14.cis.out,a14.trans.out,a8.cis.out,a8.trans.out,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/meta.cistrans.subsets.RData")
#}
#all<-rbind(a8.cis.out,a14.cis.out,a8.trans.out,a14.trans.out)
#library(ggplot2)

#p1<-ggplot(all, aes(x=HetISq)) + 
#geom_density(aes(group=subset,colour=subset))
#ggsave(p1,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/heterogeneity_meta.png",height=6,width=8)

#######

library(ggplot2)
path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis"

load(paste0(path,"/results/16/16_clumped.rdata"))
a<-clumped

f<-function(x){
y<-strsplit(x,split="")
length(which(y[[1]]=="?"))}

dir<-as.list(a$Direction)
NoDirection<-sapply(dir,f)

f<-function(x){
y<-strsplit(x,split="")
length(which(y[[1]]=="+"))}

dir<-as.list(a$Direction)
Directionplus<-sapply(dir,f)

f<-function(x){
y<-strsplit(x,split="")
length(which(y[[1]]=="-"))}

dir<-as.list(a$Direction)
Directionminus<-sapply(dir,f)

a<-data.frame(a,NoDirection,Directionplus,Directionminus)
o<-order(a$NoDirection,decreasing=T)
a<-a[o,]

a14.cis.out<-data.frame(a[which(a$cis==T & a$pval<1e-14),],subset = "cis p<1e-14")
a14.trans.out<-data.frame(a[which(a$cis==F & a$pval<1e-14),],subset = "trans p<1e-14")
a8.cis.out<-data.frame(a[which(a$cis==T & a$pval<5e-8),],subset = "cis p<5e-8")
a8.trans.out<-data.frame(a[which(a$cis==F & a$pval<5e-8),],subset = "trans p<5e-8")
all<-rbind(a8.cis.out,a14.cis.out,a8.trans.out,a14.trans.out)

p1<-ggplot(all, aes(x=HetISq)) + 
geom_density(aes(group=subset,colour=subset))
ggsave(p1,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/heterogeneity_metaclumped.png",height=6,width=8)

p1<-ggplot(a14.cis.out, aes(x=HetISq, y=Effect)) +
    geom_point() +    
    geom_smooth(method=lm,   # Add linear regression line
                se=TRUE)

ggsave(p1,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/cis1e14hetvsbeta.png",height=6,width=8)


y<-split(a,as.factor(a$Direction))
m<-sapply(y,function(x) mean(x$HetISq))
a14.cis.out$Direction<-as.factor(a14.cis.out$Direction)
a14.cis.out$HetISq<-as.numeric(a14.cis.out$HetISq)

labs <- data.frame(table(a14.cis.out$Direction))
m<-match(a14.cis.out$Direction,labs[,1])
a14.cis.out<-data.frame(a14.cis.out,Ncatdir=as.character(labs[m,-1]))


a14.cis.out$Direction<-as.factor(a14.cis.out$Direction)
a14.cis.out$Direction <- factor(a14.cis.out$Direction, levels = a14.cis.out$Direction[order(as.numeric(as.character(a14.cis.out$Ncatdir)),decreasing=T)])

p1<-ggplot(data=a14.cis.out, aes(x=a14.cis.out$Direction, y=a14.cis.out$HetISq)) +
geom_boxplot() +
theme(axis.text.x = element_text(size=8, angle=90,hjust=-0.2,vjust=0.2)) +
geom_text(data=a14.cis.out, aes(x=a14.cis.out$Direction,y = 105,label = Ncatdir),vjust = 0,size=2,angle=90) +
ylab("I2") 
ggsave(p1,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/cis1e14hetvsdir.png",height=6,width=16)

####
a14.cis.out$i2cat<-cut(as.numeric(as.character(a14.cis.out$HetISq)), breaks = seq(0,100,by=10))
w<-which(is.na(a14.cis.out$i2cat))
table(a14.cis.out$HetISq[w])
a14.cis.out$i2cat[w]<-"(0,10]"

labs <- data.frame(table(a14.cis.out$i2cat))
m<-match(a14.cis.out$i2cat,labs[,1])
a14.cis.out<-data.frame(a14.cis.out,Ncati2=as.character(labs[m,-1]))

m <- data.frame(c(by(abs(a14.cis.out$Effect), a14.cis.out$i2cat, max)))
m2<-match(a14.cis.out$i2cat,row.names(m))
a14.cis.out<-data.frame(a14.cis.out,maxbeta=abs(m[m2,]))

a14.cis.out$Effect_abs<-abs(a14.cis.out$Effect)
p1<-ggplot(a14.cis.out, aes(x=i2cat, y=Effect_abs)) +
geom_boxplot() +
geom_text(data=a14.cis.out, aes(x=a14.cis.out$i2cat,y = (maxbeta+0.1),label = Ncati2),vjust = 0,size=3) +
ylab("BETA") 
ggsave(p1,file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/cis1e14hetcatvsbeta.png",height=6,width=16)













