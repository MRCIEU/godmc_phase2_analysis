path="../results/enrichments"
l<-list.files(path=path,pattern="controls")

cis<-grep("controls_cis",l)
trans<-grep("controls_trans",l)
l1<-l[-c(cis,trans)]

cis<-l[cis]
trans<-l[trans]

controlset<-NULL
for (i in 1:10){

load(paste(path,"/controls",i,".rdata",sep=""))
cat(i,"\n")

#first remove rows with NAs only

controls2<-controls[,-1]

ind <- rowSums(is.na(controls2)) == length(controls2) 
out <- rowSums(controls2, na.rm = TRUE) 
out[ind] <- NA 

length(which(out==1))
#[1] 57978
length(which(out==0))
#[1] 57978

set<-controls$SNP[which(out==1)]
controlset<-unique(c(controlset,set))
}

