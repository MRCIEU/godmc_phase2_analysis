arguments<-commandArgs(T)
i<-as.numeric(arguments[1])
 
###
load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection.rdata")
length(unique(f.all$SNP))

#[1] 10085072
 
table(f.all$mQTL)

  FALSE    TRUE 
9885229  199843 
 
#  FALSE    TRUE 
#9885229  199843 
 
 
#path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16"
#l<-list.files(path,pattern=".txt.gz")
 
#df.out<-data.frame()
#for (i in 1:length(l)){
#cat(i,"\n")
#m<-read.table(paste(path,l[i],sep="/"),he=T)
#spl<-strsplit(as.character(m$MarkerName),split="_")
#spl<-do.call("rbind",spl)
#df<-data.frame(spl,m$Pvalue)
#df.out<-rbind(df.out,df)
#}
#names(df.out)<-c("snp","cpg","pval")
#save(file="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/meta_pval.Robj")
#load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/meta_pval.Robj")
 
#dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16"
#cd $dir
#touch snpcpgpval.txt
#touch tmp.txt
#for i in `seq 1 962`;do
#echo $i
#zcat /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_$i.txt.gz |awk 'NR>1 {print $1,$8, '$i'}' >>snpcpgpval.txt
#zcat /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_$i.txt.gz |awk -F'[ :_]' 'NR>1 {print $1,$2,$3}' >>tmp.txt
#done
 
#paste -d ' ' tmp.txt snpcpgpval.txt > snpcpgpval.tmp
#mv snpcpgpval.tmp /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.txt
#rm /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/tmp.txt
#gzip /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.txt
 
#for i in `seq 1 23`; do
#echo $i
#zcat $dir/snpcpgpval.txt.gz| awk -v mychr="chr$i" '$1==mychr {print $0}' > $dir/snpcpgpval.chr$i.txt
#done
 
df.out<-read.table(paste("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".txt",sep=""),he=F)
spl<-strsplit(as.character(df.out$V4),split="_")
spl<-do.call("rbind",spl)
df.out<-data.frame(cpg=spl[,2],df.out)
 
retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")

#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
  
df.out<-df.out[which(df.out$cpg%in%retaincpg),]

library(meffil)

y<-meffil.get.features("450k")
m<-match(df.out$cpg,y$name)
 
df.out<-data.frame(df.out,y[m,c("chromosome","position")])
names(df.out)<-c("cpg","snpchr","snppos","type","id","pval","chunk","cpgchr","cpgpos")
w<-which(df.out$snpchr==df.out$cpgchr&abs(df.out$snppos-df.out$cpgpos)<100000)
df.out$cis<-"FALSE"
df.out$cis[w]<-"TRUE"
df.out$trans<-"FALSE"
w<-which(df.out$snpchr==df.out$cpgchr&abs(df.out$snppos-df.out$cpgpos)>100000)
df.out$trans[w]<-"TRUE"
w<-which(df.out$snpchr!=df.out$cpgchr)
df.out$trans[w]<-"TRUE"

#cl[which(cl$snp=="chr9:98252950:SNP"),]
 
#ALL-cis+trans
o<-order(df.out$pval)
df.out2<-df.out[o,]
 
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
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/mqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/mqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/mqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/mqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,sep="")

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/annotation_selection/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)

##
#trans
#set cispval on 1
w<-which(df.out$cis==TRUE & df.out$trans==FALSE)
df.out$pval[w]<-1

o<-order(df.out$pval)
df.out2<-df.out[o,]
 
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
f.all2<-f.all.chr[m,c("snppos","ihs","fst","xpehhchb","xpehhyri")]

w<-which(is.na(f.all2$ihs))
f.all2$ihs[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/transmqtl_ihs/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$fst))
f.all2$fst[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/transmqtl_fst/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhyri))
f.all2$xpehhyri[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/transmqtl_xpehhchb/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

w<-which(is.na(f.all2$xpehhchb))
f.all2$xpehhchb[w]<-0
df.out3<-df.out2[-w,]
write.table(df.out3,paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/pval/transmqtl_xpehhyri/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2$annot<-paste(f.all2$ihs,f.all2$fst,f.all2$xpehhchb,f.all2$xpehhyri,sep="")

write.table(data.frame(f.all2$snppos,f.all2$annot),paste("/panfs/panasas01/shared-godmc/GARFIELD/garfield-data/annotation_selection_trans/chr",i,sep=""),sep=" ",col.names=F,row.names=F,quote=F)


















