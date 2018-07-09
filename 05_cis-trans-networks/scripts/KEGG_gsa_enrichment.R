library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stats)
library(missMethyl)
library(KEGGREST)

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
go.res[[i]] <- gometh(sig.cpg=sigCpGs, all.cpg=unique(communities$cpg),collection="KEGG")
print(table(go.res[[i]]$FDR<0.05))
}
}
save(go.res,file="../results/KEGGenrichments.rdata")


rm<-c(149,159)
clusters<-clusters[-rm]

l.out<-NULL
for (i in 1:length(go.res)){
l<-length(which(go.res[[i]]$FDR<0.05))
l.out<-rbind(l.out,l)
}
w<-which(l.out>0) #27 clusters have FDR<0.05 for GO enrichment

clusters<-clusters[w]

df.out<-data.frame()
for (j in 1:length(w)){
k<-w[j]
w1<-which(go.res[[k]]$FDR<0.05)
df<-data.frame(cluster=clusters[j],go.res[[k]][w1,])
df.out<-rbind(df.out,df)
}
save(df.out,file="../results/KEGGenrichments_fdr0.05.rdata")

t<-data.frame(table(df.out$cluster))
t<-t[t$Freq>0,]

for (i in 1:nrow(t)) {
test<-df.out[which(df.out$cluster==t$Var1[i]),]
o<-order(test$FDR)
test<-test[o,]
print(test[1:5,])
}

#7: PI3K-Akt signaling pathway
#8: SLE, Alcoholism, Viral carcinogenesis
#9: Alzheimer's disease
#11: Chemokine signaling pathway
#19: Phosphatidylinositol signaling system
#21: SLE
#23: Amino sugar and nucleotide sugar metabolism
#33: Focal adhesion, insulin signaling pathway
#36: TNF signaling pathway, PI3K-Akt signaling pathway
#39: Th1 and Th2 cell differentiation
#61: Insulin secretion, PI3K-Akt signaling pathway
#67: Antigen processing and presentation
#74: Neurotrophin signaling pathway
#109: RNA degradation
#110: Glutathione metabolism
#122: Antigen processing and presentation
#162: Aldosterone-regulated sodium reabsorption
#182: Ras signaling pathway
#219: Kaposi's sarcoma-associated herpesvirus infection
#220: cell cycle
#227: Glycolysis / Gluconeogenesis
#314: Type II diabetes mellitus 
df<-data.frame(table(df.out$Pathway))
df[df$Freq>3,]


df.out[grep("path:hsa04151",row.names(df.out)),]
#               cluster                    Pathway   N DE         P.DE
#path:hsa04151        7 PI3K-Akt signaling pathway 103 11 8.144342e-08
#path:hsa041512      36 PI3K-Akt signaling pathway 103  8 3.067687e-06
#path:hsa041511      33 PI3K-Akt signaling pathway 103  5 4.954816e-04
#path:hsa041513      61 PI3K-Akt signaling pathway 103  4 4.539063e-04
#                        FDR
#path:hsa04151  1.319383e-05
#path:hsa041512 1.987861e-04
#path:hsa041511 4.267595e-03
#path:hsa041513 4.902188e-02

df.out[which(df.out$Pathway=="Rap1 signaling pathway"),]
#               cluster                Pathway  N DE         P.DE         FDR
#path:hsa040152      33 Rap1 signaling pathway 79  5 0.0001862706 0.002623987
#path:hsa04015        7 Rap1 signaling pathway 79  6 0.0006607933 0.007136567
#path:hsa040151      23 Rap1 signaling pathway 79  5 0.0004659976 0.025163868
#path:hsa040154     182 Rap1 signaling pathway 79  3 0.0003706541 0.030022985
#path:hsa040153      36 Rap1 signaling pathway 79  4 0.0059795522 0.037987743

df.out[which(df.out$Pathway=="Ras signaling pathway"),]
#               cluster               Pathway  N DE         P.DE          FDR
#path:hsa040143      36 Ras signaling pathway 68  6 0.0000312465 0.0007231332
#path:hsa04014        7 Ras signaling pathway 68  5 0.0020036272 0.0154565525
#path:hsa040141      23 Ras signaling pathway 68  5 0.0002336534 0.0189259219
#path:hsa040144     182 Ras signaling pathway 68  3 0.0002047078 0.0221084459
#path:hsa040142      33 Ras signaling pathway 68  3 0.0101165153 0.0399725727

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



