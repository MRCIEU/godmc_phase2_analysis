library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stats)
library(missMethyl)

load("../results/communities.rdata")
nrow(communities) #13651
spl<-do.call("rbind",strsplit(communities$snp,split=":"))
w<-which(spl[,1]=="chr6"&spl[,2]>29570005&spl[,2]<33377657)
communities<-communities[-w,]
nrow(communities) #12405

communities<-data.frame(communities)
#clusters<-c("8","31","12")
df<-data.frame(table(communities$cluster)) #2316
clusters<-df[which(df$Freq>10),1] #165

go.res<-list()
for (i in 1:length(clusters)){
cat(i,"\n")
if((i!=149) && (i!=159)){
sigCpGs<-unique(as.character(communities[which(communities$cluster==clusters[i]),1]))
go.res[[i]] <- gometh(sig.cpg=sigCpGs, all.cpg=unique(communities$cpg))
print(table(go.res[[i]]$FDR<0.05))
}}
save(go.res,file="../results/goenrichments.rdata")

rm<-c(149,159)
clusters<-clusters[-rm]

l.out<-NULL
for (i in 1:length(go.res)){
l<-length(which(go.res[[i]]$FDR<0.05))
l.out<-rbind(l.out,l)
}
w<-which(l.out>0) #11 clusters have FDR<0.05 for GO enrichment

clusters<-clusters[w]

df.out<-data.frame()
for (j in 1:length(w)){
k<-w[j]
w1<-which(go.res[[k]]$FDR<0.05)
df<-data.frame(cluster=clusters[j],go.res[[k]][w1,])
df.out<-rbind(df.out,df)
}
save(df.out,file="../results/goenrichments_fdr0.05.rdata")

t<-table(df.out$cluster)
t[t>0]

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

####
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(org.Hs.eg.db)
library(limma)
library(missMethyl)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

sigcpgs<-communities[communities$cluster==8,"cpg"]
allcpgs<-unique(communities$cpg)

mappedEz <- getMappedEntrezIDs(sigcpgs,allcpgs,array.type="450k")
mappedEz$sig.eg


#Community 149 and 159 give the following error.
#Error in BiasedUrn::pWNCHypergeo(S[i, 1 + j], S[i, "N"], Total - S[i,  : 
#  Invalid value for odds

#g8[g8$FDR<0.05,]
#                                               Term Ont  N DE         P.DE
#GO:0006323                            DNA packaging  BP 49  5 2.381259e-05
#GO:0006333        chromatin assembly or disassembly  BP 47  5 1.295659e-05
#GO:0006334                      nucleosome assembly  BP 37  5 4.168682e-06
#GO:0031497                       chromatin assembly  BP 40  5 6.168938e-06
#GO:0034728                  nucleosome organization  BP 44  5 1.352560e-05
#GO:0065004             protein-DNA complex assembly  BP 50  5 1.137304e-05
#GO:0071824 protein-DNA complex subunit organization  BP 57  5 2.995467e-05
#GO:0000786                               nucleosome  CC 30  5 5.566210e-07
#GO:0032993                      protein-DNA complex  CC 49  6 1.145897e-06
#GO:0044815                    DNA packaging complex  CC 31  5 5.979691e-07
#                   FDR
#GO:0006323 0.037274634
#GO:0006333 0.023818578
#GO:0006334 0.014682099
#GO:0031497 0.017381600
#GO:0034728 0.023818578
#GO:0065004 0.023818578
#GO:0071824 0.042200138
#GO:0000786 0.004212094
#GO:0032993 0.005381133
#GO:0044815 0.004212094

#[1] Term Ont  N    DE   P.DE FDR 
#<0 rows> (or 0-length row.names)
#[1] Term Ont  N    DE   P.DE FDR 
#<0 rows> (or 0-length row.names)

