arguments<-commandArgs(T)
i<-as.numeric(arguments[1])


#load("../results/enrichments/snpcontrolsets.rdata")

#r<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/iHS_CEU.whole_genome.pvalues.gz",he=T)
#names(r)[4:5]<-c("iHS_score","iHS_pval")
#r2<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/FstGLOB_CEU_u_YRI_u_CHB.whole_genome.pvalues.gz",he=T)
#names(r2)[4:5]<-c("Fst_score","fst_pval")

#m<-match(f.all$SNP,r$snpID)
#m2<-match(f.all$SNP,r2$snpID)
#f.all<-data.frame(f.all,r[m,c("iHS_score","iHS_pval")],r2[m2,c("Fst_score","fst_pval")])

#save(f.all,file="../results/enrichments/snpcontrolsets_selection.rdata")

#test<-f.all[,c("snpchr","snppos","iHS_pval")]
#w<-which(f.all.chr$iHS_pval<0.05)
#w<-which(test$iHS_pval<0.05)
#length(w)
#[1] 811651
#test$iHS_pval[w]<-0
#w<-which(test$iHS_pval>0.05)
#test$iHS_pval[w]<-1
#table(test$snpchr,test$iHS_pval)
       
 #            0      1
 # chr1   63129 503117
 # chr2   67366 545807
 # chr3   58037 468467
 # chr4   58329 479556
 # chr5   51689 415769
 # chr6   52660 420578
 # chr7   46826 379099
 # chr8   44036 364423
 # chr9   36323 284185
 # chr10  41612 331184
 # chr11  40003 316514
 # chr12  38965 311456
 # chr13  31428 246479
 # chr14  26616 214017
 # chr15  23038 180298
 # chr16  26185 201739
 # chr17  22390 169127
 # chr18  23903 191085
 # chr19  16834 129773
 # chr20  18443 147753
 # chr21  11851  91957
 # chr22  11988  87044
 # chr23      0      0


###
load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/enrichments/snpcontrolsets_selection.rdata")
length(unique(f.all$SNP))
#[1] 10085072

table(f.all$mQTL)

#  FALSE    TRUE 
#9885229  199843 
w1<-which(f.all$iHS_pval>0.05)
w2<-which(f.all$iHS_pval<0.05)
f.all$iHS_pval[w2]<-1
f.all$iHS_pval[w1]<-0

w1<-which(f.all$fst_pval>0.05)
w2<-which(f.all$fst_pval<0.05)
f.all$fst_pval[w2]<-1
f.all$fst_pval[w1]<-0

table(f.all$snpchr,f.all$iHS_pval)
table(f.all$snpchr,f.all$fst_pval)

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

o<-order(df.out$pval)
df.out<-df.out[o,]

bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.bim")
bim<-bim[bim$V1==i,]
m<-match(bim$V2,df.out$snppos)
df.out2<-data.frame(snp=bim$V4,pval=df.out$snppos[m])
w<-which(is.na(df.out2$pval))
df.out2$pval[w]<-1

length(which(!is.na(df.out2$pval)))


chr<-paste("chr",i,sep="")
f.all.chr<-f.all[f.all$snpchr==chr,]
m<-match(df.out2$snp,f.all.chr$snppos)
f.all2<-f.all.chr[m,c("snppos","iHS_pval","fst_pval")]

w<-which(is.na(f.all2$iHS_pval))
df.out2<-df.out2[-w,]
write.table(df.out2,paste("/panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/pval/mqtl/chr",i,sep=""),sep=" ",quote=F,col.names=F,row.names=F)

f.all2<-f.all2[-w,]
write.table(f.all2[,1:2],paste("~/GARFIELD/garfield-data/sel_annotation/chr",j,sep=""),sep=" ",col.names=F,row.names=F,quote=F)




#/panfs/panasas01/sscm/epzjlm/GARFIELD/garfield
#
#for i in `seq 1 23`; do
#echo $i
#./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl/chr$i ../garfield-data/annotation/chr$i -1 > ../garfield-data/prep/chr$i
#done