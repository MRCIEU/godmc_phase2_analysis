path="/projects/MRC-IEU/groups/godmc/sftp/GoDMC/*/"
l1<-list.dirs("/panfs/panasas01/shared-godmc/sftp/GoDMC")
g<-grep("/results/01",l1)
l1<-l1[g]
g<-grep("logs",l1)
l1<-l1[-g]
g<-grep("mbustamante/results/01",l1)
l1<-l1[-g]

spl<-strsplit(l1,split="/")
spl<-do.call("rbind",spl)
#l1<-l1[-1]
#l2<-list.files(l,pattern="_01.tgz")

for (i in 1:length(l1)){
load(paste(l1[i],"methylation_summary.RData",sep="/"))
names(meth_summary)<-paste(names(meth_summary),sep=".")

if(length(l)>0){

for (j in 1:length(l)){
r<-read.table(paste(l1[i],l[j],sep="/"))

}

}
}