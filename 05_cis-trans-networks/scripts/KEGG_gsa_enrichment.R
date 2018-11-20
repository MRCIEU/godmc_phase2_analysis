library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stats)
library(missMethyl)
library(KEGGREST)
library(meffil)

#load("../results/communities.rdata")
#nrow(communities) #13651
#spl<-do.call("rbind",strsplit(communities$snp,split=":"))
#w<-which(spl[,1]=="chr6"&spl[,2]>24570005&spl[,2]<38377657)
#w<-which(spl[,1]=="chr6"&spl[,2]>24570005&spl[,2]<38377657)

#communities<-communities[-w,]
#nrow(communities) #11898

#communities<-data.frame(communities)
#clusters<-c("8","31","12")
#df<-data.frame(table(communities$cluster)) #2316
#clusters<-df[which(df$Freq>10),1] #154
#rm<-c(138,148)
#clusters<-clusters[-rm]

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
t<-data.frame(table(mem$cluster))
t[t$Freq>50,]
#  Var1 Freq
#1    1   53
#2    2  253
#6    6   64

y<-meffil.get.features()
y<-y[,c("name","chromosome","position")]
m<-match(mem$cpg,y$name)
mem<-data.frame(mem,y[m,])

hla<-which(mem$chromosome=="chr6"&mem$position>2470005&mem$position<38377657)
length(hla) #111
mem2<-mem[-hla,]

t2<-data.frame(table(mem2$cluster))
m<-merge(t,t2,by.x=1,by.y=1)
m<-m[order(as.numeric(m[,1])),]

clusters<-t2[which(t2$Freq>=10),1] #55

go.res<-list()
for (i in 1:length(clusters)){
cat(i,"\n")
cat(clusters[i],"\n")
sigCpGs<-unique(as.character(mem2[which(mem2$cluster==clusters[i]),1]))
go.res[[i]] <- gometh(sig.cpg=sigCpGs, all.cpg=unique(mem2$cpg),collection="KEGG")
#sigCpGs<-unique(as.character(communities[which(communities$cluster==clusters[i]),1]))
#go.res[[i]] <- gometh(sig.cpg=sigCpGs, all.cpg=unique(communities$cpg),collection="KEGG")
go.res[[i]]$cluster<-clusters[i]
go.res[[i]]$KEGG<-row.names(go.res[[i]])
print(table(go.res[[i]]$FDR<0.05))
}

save(go.res,file="../results/KEGGenrichments.rdata")

df.out<-do.call("rbind",go.res)
df.out<-df.out[which(df.out$FDR<0.05),]

y<-.Machine$double.xmin #2.225074e-308
w<-which(df.out$FDR==0)
df.out$FDR[w]<-y
w<-which(df.out$P.DE==0)
df.out$P.DE[w]<-y

save(df.out,file="../results/KEGGenrichments_fdr0.05.rdata")
write.table(df.out,"kegg_enrichments_fdr0.05.txt",sep="\t",quote=F,row.names=F,col.names=T)

df.out[which(df.out$FDR<1e-300),]
t<-data.frame(table(df.out$cluster))
t<-t[t$Freq>0,]

#   Var1 Freq
#2     2  110
#3     3    2
#5     5    8
#9     9    1
#18   18    3
#19   19    7
#21   21    1
#22   22   33
#25   25    2
#26   26    1
#64   64    1


for (i in 1:nrow(t)) {
test<-df.out[which(df.out$cluster==t$Var1[i]),]
o<-order(test$FDR)
test<-test[o,]
print(test[1:5,])
}


df<-data.frame(table(df.out$Pathway))
df[df$Freq>2,]

df.out[which(df.out$Pathway=="Acute myeloid leukemia"),]
#                               Pathway  N DE         P.DE        FDR cluster
#path:hsa052211  Acute myeloid leukemia 15  3 0.0016727279 0.00746999       2
#path:hsa0522117 Acute myeloid leukemia 15  2 0.0004382108 0.04763808      18
#path:hsa0522121 Acute myeloid leukemia 15  2 0.0014552795 0.02499668      22
#                         KEGG
#path:hsa052211  path:hsa05221
#path:hsa0522117 path:hsa05221
#path:hsa0522121 path:hsa05221

df.out[which(df.out$Pathway=="Ras signaling pathway"),]
#                              Pathway  N DE         P.DE          FDR cluster
#path:hsa040141  Ras signaling pathway 45  6 6.503671e-05 0.0005889435       2
#path:hsa040142  Ras signaling pathway 45  3 2.680502e-04 0.0436921904       3
#path:hsa0401420 Ras signaling pathway 45  3 7.584069e-05 0.0247240633      21
#                         KEGG
#path:hsa040141  path:hsa04014
#path:hsa040142  path:hsa04014
#path:hsa0401420 path:hsa04014




#2: transcription
#8: chromatin assembly
#20 embryonic skeletal system development
#25: transcription
#41: norepinephrine secretion
#86: immune response
#122: immune response
#126: circadian behavior
#202: embryonic morphogenesis
#258: thyroid hormone generation
#298: cell adhesion


#     V1          V2
#6     8 0.013888889
#11   14 0.023255814
#26   32 0.071428571
#29   36 0.007633588
#52   68 0.052631579
#70  105 0.090909091
#86  134 0.030303030
#98  162 0.125000000
#103 182 0.013157895
#118 214 0.043478261
#119 215 0.068181818

