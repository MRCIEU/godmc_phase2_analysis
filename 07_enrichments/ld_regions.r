library(tidyverse)
library(dplyr)
library(meffil)
library(GenomicRanges)

# Get clumped data
load("../results/16/16_clumped.rdata")
cl<-data.frame(clumped)
cl.out<-data.frame()

y<-meffil.get.features("450k")

for (i in 1:23){
cat(i,"\n")
cl_chr<-cl[which(cl$snpchr%in%paste0("chr",i)),]
a <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.",i,".ld.formatted.gz"))
a<-data.frame(a)

agg1<-a %>%  
  group_by(BP_A) %>%  
  dplyr::summarize (bpa = first(BP_A),start_bp = min(BP_B),  
             stop_bp = max(BP_B),nproxies = n(),chr=first(CHR_A))

agg2<-a %>%  
  group_by(BP_B) %>%  
  dplyr::summarize (bpb = first(BP_B),start_bp = min(BP_A),  
             stop_bp = max(BP_A),nproxies = n(),chr=first(CHR_B))

b <- full_join(agg1, agg2, c("BP_A" = "BP_B"))

#b$final_min <- b$start_bp.x
#b$final_max <- b$stop_bp.y
#b$final_min[is.na(b$final_min)] <- b$bpa[is.na(b$final_min)]
#b$final_max[is.na(b$final_max)] <- b$bpb[is.na(b$final_max)]
#b<-data.frame(b)

b2<-data.frame(b$nproxies.x,b$nproxies.y)
b$nproxies<-apply(b2,1,function(x) y<-sum(x,na.rm=T))

b2<-data.frame(b$start_bp.x,b$start_bp.y)
b$final_min<-apply(b2,1,function(x) y<-min(x,na.rm=T))

b2<-data.frame(b$stop_bp.x,b$stop_bp.y)
b$final_max<-apply(b2,1,function(x) y<-max(x,na.rm=T))


b$CHR<-b$chr.x
b$CHR[is.na(b$CHR)]<-b$chr.y[is.na(b$CHR)]
b$CHR<-paste0("chr",b$CHR)

table(is.na(b$final_min))
table(is.na(b$final_max))
table(b$final_min <= b$final_max)
table(b$final_min <= b$final_max)


b <- dplyr::select(b, CHR=CHR,SNP=BP_A, min=final_min, max=final_max,nproxies=nproxies)
summary(b$max - b$min)

#add MAF

f <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/frq/eur.filtered.",i,".formatted.gz"))
spl<-strsplit(f$SNP,split=":")
spl<-do.call("rbind",spl)
colnames(spl)<-c("snpchr","snppos","snptype")
f<-data.frame(spl,f)
m<-match(f$snppos,b$SNP)
f<-data.frame(snpchr=f[,c("snpchr")],b[m,-1:-2],f[,c("snppos","MAF","snptype","SNP")])
f$snppos<-as.numeric(as.character(f$snppos))
f$min<-as.numeric(as.character(f$min))
f$max<-as.numeric(as.character(f$max))

f$min[is.na(f$min)]<-f$snppos[is.na(f$min)]
f$max[is.na(f$max)]<-f$snppos[is.na(f$max)]
f$nproxies[is.na(f$nproxies)]<-0

#add distance to CpG

#load("../03_clumping_16/cpg_pos.rdata")
#cpgpos<-cpgpos[which(cpgpos$cpgchr==paste0("chr",i)),]

cpgpos<-y[which(y$chromosome==paste0("chr",i)),]
cpgpos$position<-as.numeric(as.character(cpgpos$position))
snppos<-as.numeric(as.character(f$snppos))

#f1<-function(x){
#x<-f$snppos
#m<-min(abs(x-cpgpos$position),na.rm=T)
#w1<-which(cpgpos$position==x+m)
#w2<-which(cpgpos$position==x-m)
#w<-na.omit(c(w1,w2))
#return(cpgpos[w,"name"])
#}

#out<-apply(f,1,f1)

closestcpg<-ddply(f, .(snppos), summarise,
        closestcpgposmin = snppos-min(abs(snppos-cpgpos$position),na.rm=T),closestcpgposmax = snppos+min(abs(snppos-cpgpos$position),na.rm=T) )

m1<-match(closestcpg[,2],cpgpos$position)
m2<-match(closestcpg[,3],cpgpos$position)
cpgpos1<-data.frame(cpgpos[m1,"name"],cpgpos[m2,"name"])

cpgpos2<-as.character(cpgpos1[,1])
#cpgpos2[!is.na(cpgpos1[,1])]<-as.character(cpgpos1[!is.na(cpgpos1[,1]),1])
cpgpos2[!is.na(cpgpos1[,2])]<-as.character(cpgpos1[!is.na(cpgpos1[,2]),2])

m<-match(cpgpos2,cpgpos$name)
dist<-abs(cpgpos[m,c("position")]-f$snppos)
f<-data.frame(f,closest450kcpg=cpgpos[m,c("name")],closest450kcpgpos=cpgpos[m,c("position")],closest450kdistance=dist) #16050654      72612

cl_chr.cis<-cl_chr[which(cl_chr$cis==T),]
cl_chr.trans<-cl_chr[which(cl_chr$cis==F),]

m<-which(f$snppos%in%cl_chr.cis$snppos)
f$cismQTL<-"FALSE"
f$cismQTL[m]<-"TRUE"

m<-which(f$snppos%in%cl_chr.trans$snppos)
f$transmQTL<-"FALSE"
f$transmQTL[m]<-"TRUE"

m<-which(f$snppos%in%cl_chr$snppos)
f$mQTL<-"FALSE"
f$mQTL[m]<-"TRUE"

mafcat<-cut(f$MAF,breaks=seq(0,0.5,0.05))
distcat<-cut(f$closest450kdistance,breaks=seq(0,10000000,10000))

w<-which(f$closest450kdistance==0)
distcat[w]<-levels(distcat)[1]

proxycat<-cut(f$nproxies,breaks=seq(0,100000,10))
w<-which(f$nproxies==0)
proxycat[w]<-levels(proxycat)[1]

f<-data.frame(f,mafcat,distcat,proxycat)

f$groups<-paste(f$mafcat,f$distcat,f$proxycat)

#check<-table(f$groups,f$mQTL)
#check<-check[check[,2]!=0,]
#check<-data.frame(group=row.names(check),nomQTL=check[,1],mQTL=check[,2])
#o<-order(check[,1])
#check[o,]

library(ggplot2)
#p1<-ggplot(x=groups,data=f, geom="bar", stat="identity",position="dodge")

labs <- data.frame(table(f$groups))
m<-match(f$groups,labs[,1])
f<-data.frame(f,Ncatgroups=as.character(labs[m,-1]))

f$groups<-as.factor(f$groups)
f$groups <- factor(f$groups, levels = f$groups[order(as.numeric(as.character(f$Ncatgroups)),decreasing=T)])


#select groups with mQTLs
t<-table(f$groups,f$mQTL)
t<-t[t[,2]>0,]
f2<-f[which(f$groups%in%row.names(t)),]

p1<-ggplot(f2, aes(x=groups),color=mQTL,fill=mQTL) + 
geom_bar(position="dodge")
ggsave(p1,file="groupsizes.pdf",width=4,height=4)

p1<-ggplot(f2, aes(x=groups),color=mQTL,fill=mQTL) + 
geom_bar(position="dodge") +
ylim(0,1000)
ggsave(p1,file="groupsizes1000.pdf",width=4,height=4)



pdf("groupsizes.pdf",height=6,width=6)
par(mfrow=c(2,2))
hist(check[,1],main=NULL,xlab="no mQTL",breaks)
hist(check[,2],main=NULL,xlab="mQTL")
dev.off()


#a2<-a[which(a$SNP_A%in%cl_chr$snp),]

#test<-a3[a3$BP_A=="16855618",]
#min(test$BP_B)
#[1] 16855677
#max(test$BP_B)
#[1] 16866339

#w<-which(agg$BP_A%in%cl_chr$snppos)
#> length(w)
#[1] 3275

#length(which(is.na(cl_chr$BP_A)))
#[1] 1208

#m<-match(cl_chr$snppos,f$snppos)
#cl_chr<-data.frame(cl_chr,snppos[m,])

cl.out<-rbind(cl.out,cl_chr)

}

w<-which(is.na(cl.out$start_bp))
cl.out$start_bp[w]<-cl.out$snppos[w]
cl.out$stop_bp[w]<-cl.out$snppos[w]

w<-which(cl.out$start_bp>cl.out$snppos)
cl.out$start_bp[w]<-cl.out$snppos[w]

w<-which(cl.out$stop_bp<cl.out$snppos)
cl.out$stop_bp[w]<-cl.out$snppos[w]

save(cl.out,file="../results/enrichments/clumpedwithldregion.rdata")


