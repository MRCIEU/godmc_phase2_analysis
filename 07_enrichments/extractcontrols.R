path="../results/enrichments"
l<-list.files(path=path,pattern="controlslist")

cis<-grep("controlslist_cis",l)
trans<-grep("controlslist_trans",l)
l1<-l[-c(cis,trans)]

cis<-l[cis]
trans<-l[trans]

load("../results/enrichments/snpcontrolsets.rdata")
####

controlset<-list()
mqtlset<-NULL

for (i in 1:10){

load(paste(path,l1[i],sep="/"))
mqtlset<-unique(c(mqtlset,unlist(mqtllist)))

snp<-unlist(controllist)
controlset[[i]]<-f.all[which(f.all$SNP%in%snp),]
table(controlset[[1]][,"mQTL"])
cat(i,"\n")
}

mqtlset<-f.all[which(f.all$SNP%in%mqtlset),]

dim(mqtlset)
#[1] 232670
dim(controlset[[1]])

save(mqtlset,controlset,file="snpsetsforLOLA.RData")

###

controlset<-list()
mqtlset<-NULL

for (i in 1:10){

load(paste(path,cis[i],sep="/"))
mqtlset<-unique(c(mqtlset,unlist(mqtllist)))

snp<-unlist(controllist)
controlset[[i]]<-f.all[which(f.all$SNP%in%snp),]
table(controlset[[1]][,"mQTL"])
table(controlset[[1]][,"cismQTL"])
table(controlset[[1]][,"transmQTL"])
cat(i,"\n")
}

mqtlset<-f.all[which(f.all$SNP%in%mqtlset),]

dim(mqtlset)
#[1] 232670
dim(controlset[[1]])

save(mqtlset,controlset,file="snpsetsforLOLAcis.RData")

##

controlset<-list()
mqtlset<-NULL

for (i in 1:10){

load(paste(path,trans[i],sep="/"))
mqtlset<-unique(c(mqtlset,unlist(mqtllist)))

snp<-unlist(controllist)
controlset[[i]]<-f.all[which(f.all$SNP%in%snp),]
table(controlset[[1]][,"mQTL"])
table(controlset[[1]][,"cismQTL"])
table(controlset[[1]][,"transmQTL"])
cat(i,"\n")
}

mqtlset<-f.all[which(f.all$SNP%in%mqtlset),]

dim(mqtlset)
#[1] 232670
dim(controlset[[1]])

save(mqtlset,controlset,file="snpsetsforLOLAtrans.RData")












