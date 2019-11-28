res.out<-data.frame()
r<-read.csv("CpGlist_CARMELI.csv")

for (i in 1:962){
cat(i,"\n")
p<-paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",i,".rdata")
load(p)
res<-res[which(res$cpg%in%r$ID),]
res.out<-rbind(res.out,res)
}
write.table(res.out,"probes_36cohorts16_unclumped.txt",sep="\t",col.names=T,row.names=F,quote=F)

length(unique(res.out$cpg))
#10902

i<-intersect(unique(res.out$cpg),r$ID)
length(i)