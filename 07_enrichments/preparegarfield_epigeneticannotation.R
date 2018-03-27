arguments<-commandArgs(T)
i<-as.numeric(arguments[1])

df.out<-read.table(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans.txt.gz",sep=""),he=T)

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")

#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(bim$V2%in%flip[,1])
#bim<-bim[-w,]

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(bim$V2%in%indels[,1])
#bim<-bim[-w,]

bim<-bim[bim$V1==i,]
m<-match(bim$V4,df.out$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out$min_pval[m])

w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1

length(which(!is.na(df.out2$pval)))

length(which(df.out2$pval<1e-14))
#[1] 65007

o<-order(as.numeric(as.character(df.out2$snp)))
df.out2<-df.out2[o,]

write.table(df.out2,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_epigenetic/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)


#TRANS ONLY -set cis only to pval of 1
w<-which(df.out$snp_cis!="FALSE")
#table(df.out[w,"cis"],df.out[w,"trans"])
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

o<-order(as.numeric(as.character(df.out2$snp)))
df.out2<-df.out2[o,]

write.table(df.out2,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_trans_epigenetic/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

#cis only
w<-which(df.out$snp_cis!="TRUE")
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

o<-order(as.numeric(as.character(df.out2$snp)))
df.out2<-df.out2[o,]

write.table(df.out2,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_cis_epigenetic/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

#cis+trans only
w<-which(df.out$snp_cis!="ambivalent")
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

o<-order(as.numeric(as.character(df.out2$snp)))
df.out2<-df.out2[o,]

write.table(df.out2,paste("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval/mqtl_ambivalent_epigenetic/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

#/panfs/panasas01/sscm/epzjlm/GARFIELD/garfield
#
#for i in `seq 1 22`; do
#echo $i
#./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl/chr$i ../garfield-data/annotation/chr$i -1 > ../garfield-data/prep/chr$i
#done

#cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep
#touch mqtlinputfile.txt

#for i in `seq 1 22`; do
#echo $i
#cat chr$i >> mqtlinputfile.txt
#done


#garfield-perm -n 10000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep/ \
#-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_iHS -g -m -t 0.0001 -progress



