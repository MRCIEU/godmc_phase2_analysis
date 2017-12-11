arguments<-commandArgs(T)
i<-as.numeric(arguments[1])
 
###
load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection.rdata")
length(unique(f.all$SNP))

#[1] 10085072
 
table(f.all$mQTL)
 
#  FALSE    TRUE 
#9852470  232602 

 
df.out<-read.table(paste("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans.txt.gz",sep=""),he=T)

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")

flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
w<-which(bim$V2%in%flip[,1])
bim<-bim[-w,]

indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
w<-which(bim$V2%in%indels[,1])
bim<-bim[-w,]

 #ALL-cis+trans
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out$min_pval[m])

w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri","sds")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$sds))
f.all2$sds[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_sds/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,f.all2$sds,sep="")

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))

y<-nchar(f.all2$annot)
m<-max(nchar(f.all2$annot))
m<-which(y==m)
f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_selection/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)

##
#trans
#set cispval on 1

w<-which(df.out$snp_cis!="FALSE")
table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri","sds")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$sds))
f.all2$sds[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/transmqtl_sds/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,f.all2$sds,sep="")

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))

y<-nchar(f.all2$annot)
m<-max(nchar(f.all2$annot))
m<-which(y==m)
f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_selection_trans/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)


#
#cis
#set trans and ambivalent pval on 1

w<-which(df.out$snp_cis!="TRUE")
table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri","sds")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$sds))
f.all2$sds[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/cismqtl_sds/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,f.all2$sds,sep="")

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))

y<-nchar(f.all2$annot)
m<-max(nchar(f.all2$annot))
m<-which(y==m)
f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_selection_cis/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)

###

#
#ambivalent
#set trans and cis pval on 1

w<-which(df.out$snp_cis!="ambivalent")
table(df.out[w,"cis"],df.out[w,"trans"])
df.out2<-df.out
df.out2[w,"pval"]<-1

o<-order(df.out2$pval)
df.out2<-df.out2[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out2$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out2$pval[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1
 
length(which(!is.na(df.out2$pval)))

chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri","sds")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$sds))
f.all2$sds[w]<-0
df.out3<-df.out2[-w,]

o<-order(df.out3$snp)
df.out3<-df.out3[o,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/ambivalentmqtl_sds/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,f.all2$sds,sep="")

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation/"
p<-paste("chr",i,sep="")
r<-read.table(paste(path,"chr",i,sep=""))
m<-match(r$V1,f.all2$snppos)
f.all2<-f.all2[m,]
f.all2$snppos<-r$V1
w<-which(is.na(m))

y<-nchar(f.all2$annot)
m<-max(nchar(f.all2$annot))
m<-which(y==m)
f.all2$annot[w]<-f.all2$annot[m[1]]

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation_selection_ambivalent/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)













