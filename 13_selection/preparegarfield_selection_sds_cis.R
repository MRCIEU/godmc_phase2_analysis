arguments<-commandArgs(T)
i<-as.numeric(arguments[1])
print(i)
###
load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection.rdata")
length(unique(f.all$SNP))

#[1] 10085072
 
table(f.all$mqtl_clumped)

#FALSE    TRUE 
#5980960  224648 
 
df.out<-read.table(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans.txt.gz",sep=""),he=T)
table(df.out$snp_cis)
df.out<-df.out[which(df.out$snp_cis!="FALSE"),]
nrow(df.out)
#452041

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")

#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(bim$V2%in%flip[,1])
#bim<-bim[-w,]

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(bim$V2%in%indels[,1])
#bim<-bim[-w,]

 #ALL-cis+trans
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out$min_pval[m])

w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
print(dim(df.out2))
length(which(!is.na(df.out2$pval)))

w<-which(is.na(f.all$sds))
f.all$sds[w]<-0
w<-which(is.na(f.all$xpehhchb))
f.all$xpehhchb[w]<-0
w<-which(is.na(f.all$xpehhyri))
f.all$xpehhyri[w]<-0
w<-which(is.na(f.all$fst))
f.all$fst[w]<-0
w<-which(is.na(f.all$ihs))
f.all$ihs[w]<-0

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri","sds")]
f.all2$annot<-paste(f.all2$sds,sep="")

w<-which(f.all2$sds<1)
df.out2$pval[w]<-1

print(table(f.all2$annot))
w<-which(df.out2$pval>1e-14&f.all2$annot>0)
print(length(w))
f.all2$annot[w]<-0


w<-which(is.na(f.all2$sds))
if(length(w)>0){df.out3<-df.out2[-w,]}
if(length(w)==0){df.out3<-df.out2}

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
print(dim(df.out3))

write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_sds_extremes/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

###

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]

f.all2$snppos<-r$V1
w<-which(is.na(m))

#table(f.all2$annot)

#     0      1     NA 
# 90895    702 103072 

y<-nchar(f.all2$annot)
table(y) #1 2

m<-min(nchar(f.all2$annot)) #1
m<-which(y==m)
#set missing annotations to 0
f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_sds_cis/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)



