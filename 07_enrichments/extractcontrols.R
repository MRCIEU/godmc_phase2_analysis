path="../results/enrichments"
l<-list.files(path=path,pattern="controlslist")

cis<-grep("controlslist_cis",l)
trans<-grep("controlslist_trans",l)
l1<-l[-c(cis,trans)]

cis<-l[cis]
trans<-l[trans]

controlset<-list()
mqtlset<-NULL

load("../results/enrichments/snpcontrolsets.rdata")

for (i in 1:10){

load(paste(path,l1[i],sep="/"))
mqtlset<-unique(c(mqtlset,unlist(mqtllist)))

snp<-unlist(controllist)
controlset[[i]]<-f.all[which(f.all$SNP%in%snp),]
cat(i,"\n")
}

mqtlset<-f.all[which(f.all$SNP%in%mqtlset),]

dim(mqtlset)
#[1] 186395
dim(controlset[[1]])

save(mqtlset,controlset,file="snpsetsforJohanna170724.RData")

