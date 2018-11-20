library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stats)
library(missMethyl)
library(meffil)
#load("../results/communities.rdata")
#nrow(communities) #13651
#spl<-do.call("rbind",strsplit(communities$snp,split=":"))
#w<-which(spl[,1]=="chr6"&spl[,2]>2470005&spl[,2]<38377657)
#communities<-communities[-w,]
#nrow(communities) #12405

#communities<-data.frame(communities)
#clusters<-c("8","31","12")
#df<-data.frame(table(communities$cluster)) #2316
#clusters<-df[which(df$Freq>10),1] #165 #154

#rm<-c(138,148,159)
#clusters<-clusters[-rm]

#Community 138, 148 and 159 give the following error.
#Error in BiasedUrn::pWNCHypergeo(S[i, 1 + j], S[i, "N"], Total - S[i,  : 
#  Invalid value for odds

load("../results/graph.rdata")
dim(dat)
#3573 -creg-tCpG pairs
communities <- merge(dat, mem, by.x="creg",by.y="cpg") #3573
sum(table(communities$cluster) >= 10)
#52
length(unique(communities$cluster))
#1615
cpg<-unique(c(communities$creg,communities$tcpg))
length(cpg)
#5109
head(mem)
#         cpg cluster
#1 cg05575639       3
#2 cg14872454       2
#3 cg03128025      66

dim(mem)
#5109 cpgs 
length(unique(mem$cluster))
#[1] 1615
sum(table(mem$cluster) >= 10)
#56

y<-meffil.get.features()
y<-y[,c("name","chromosome","position")]
m<-match(mem$cpg,y$name)
mem<-data.frame(mem,y[m,])

hla<-which(mem$chromosome=="chr6"&mem$position>2470005&mem$position<38377657)
length(hla) #111
mem2<-mem[-hla,]

df<-data.frame(table(mem2$cluster)) #1613
clusters<-df[which(df$Freq>=10),1] #55

go.res<-list()
for (i in 1:length(clusters)){
cat(clusters[i],"\n")
sigCpGs<-unique(as.character(mem2[which(mem2$cluster==clusters[i]),1]))
go.res[[i]] <- gometh(sig.cpg=sigCpGs, all.cpg=unique(mem2$cpg))
go.res[[i]]$cluster<-clusters[i]
go.res[[i]]$GO<-row.names(go.res[[i]])
print(table(go.res[[i]]$FDR<0.05))
}

#}
save(go.res,file="../results/goenrichments.rdata")

df.out<-do.call("rbind",go.res)
df.out<-df.out[which(df.out$FDR<0.05),]

y<-.Machine$double.xmin
w<-which(df.out$FDR==0)
df.out$FDR[w]<-y
w<-which(df.out$P.DE==0)
df.out$P.DE[w]<-y

save(df.out,file="../results/goenrichments_fdr0.05.rdata")
write.table(df.out,"go_enrichments_fdr0.05.txt",sep="\t",quote=F,row.names=F,col.names=T)