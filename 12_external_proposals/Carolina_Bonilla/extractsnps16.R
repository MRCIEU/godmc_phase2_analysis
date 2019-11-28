res.out<-data.frame()
r<-read.table("snps.txt")
for (i in 1:962){
cat(i,"\n")
p<-paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",i,".rdata")
load(p)
res<-res[which(res$snp%in%r[,1]),]
res.out<-rbind(res.out,res)
}
write.table(res.out,"snps_36cohorts16.txt",sep="\t",col.names=T,row.names=F,quote=F)
