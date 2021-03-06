path="../results/enrichments"
l<-list.files(path=path,pattern="controlslist")
g<-grep("6prop",l)
l<-l[g]

cis<-grep("controlslist6prop_cis",l)
trans<-grep("controlslist6prop_trans",l)
l1<-l[-c(cis,trans)]

cis<-l[cis]
trans<-l[trans]

load("../results/enrichments/snpcontrolsets_selection.rdata")

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
#[1] 232602
#232602
table(mqtlset$snp_cis)

#ambivalent      FALSE       TRUE 
#      2842      13516     216244 

table(controlset[[i]]$snp_cis)

#ambivalent      FALSE       TRUE 
#         0          0          0 

save(mqtlset,controlset,file="snpsetsforLOLA6prop.RData")

###

controlset<-list()
mqtlset<-NULL

for (i in 1:10){

load(paste(path,cis[i],sep="/"))
mqtlset<-unique(c(mqtlset,unlist(mqtllist)))

snp<-unlist(controllist)
controlset[[i]]<-f.all[which(f.all$SNP%in%snp),]
table(controlset[[1]][,"mQTL"])
table(controlset[[1]][,"snp_cis"])
cat(i,"\n")
}

mqtlset<-f.all[which(f.all$SNP%in%mqtlset),]

dim(mqtlset)
#[1] 232670
dim(controlset[[1]])

table(mqtlset$snp_cis)

#ambivalent      FALSE       TRUE 
#         0          0     216244 

#ambivalent      FALSE       TRUE 
#         0          0       8227 

save(mqtlset,controlset,file="snpsetsforLOLAcis6prop.RData")

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

save(mqtlset,controlset,file="snpsetsforLOLAtrans6prop.RData")












